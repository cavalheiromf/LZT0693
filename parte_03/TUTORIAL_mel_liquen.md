# Tutorial Secundário: MaAsLin3 Simplificado ("Sliced")

Este é um tutorial complementar e introdutório ao **MaAsLin3**. Antes de avançarmos para a integração complexa de múltiplos datasets (a proposta principal da Parte 03), faremos uma análise de enriquecimento simples utilizando **apenas os dados já processados da Parte 01**.

## 🎯 Objetivo
Avaliar os táxons (famílias) diferencialmente abundantes ao comparar amostras de **mel** vs **líquen**. 

Como reduzimos o escopo da análise selecionando apenas duas categorias de uma mesma corrida de sequenciamento, chamamos essa abordagem de um **"slice"** (uma "fatia" do banco de dados). O uso de um slice é uma excelente prática para:
1. Aprender a mecânica e os parâmetros básicos de uma nova ferramenta estatística de forma rápida.
2. Validar hipóteses biológicas específicas sem o ruído ou o peso computacional de amostras não relacionadas à pergunta.

---

## 1. Carregando e Preparando a "Fatia" (Slice)

Nesta etapa, você criará o script `scripts/maaslin3_sliced.R`. Note que nós não precisaremos integrar tabelas aqui; vamos ler diretamente o objeto Phyloseq que salvamos no final do módulo 01 e **filtrá-lo**.

Crie o arquivo e adicione o código abaixo:

```r
# ==============================================================================
# Script: maaslin3_sliced.R
# Objetivo: Rodar análise diferencial (mel vs liquen) usando apenas
#           um "slice" dos dados da Parte 01.
# ==============================================================================

# 1. Carregar pacotes
.libPaths(c("/opt/R/sharedLibs/4.3", .libPaths()))
library(phyloseq)
library(maaslin3)

# 2. Carregar o objeto Phyloseq da Parte 01
ps <- readRDS("../parte_01/results/phyloseq_obj.rds")
cat("Phyloseq original (Parte 01):\n")
print(ps)

# 3. Criar o "Slice": Filtrar para manter apenas Mel e Liquen
ps_sliced <- subset_samples(ps, group %in% c("mel", "liquen"))

# 4. Agrupar dados a nível de Família
ps_fam <- tax_glom(ps_sliced, taxrank = "Family", NArm = FALSE)
cat("\nPhyloseq após agrupar por Família (Apenas Mel e Liquen):\n")
print(ps_fam)

# 5. Extrair Tabela de Abundância (ASVs) e Metadados
# MaAsLin3 exige que as features (ASVs) sejam as colunas e as amostras as linhas
asv_table <- as.data.frame(t(otu_table(ps_fam)))
metadata  <- as(sample_data(ps_fam), "data.frame")

# Adicionar taxonomia aos nomes das features para o gráfico ficar legível
tax_table_df <- as.data.frame(tax_table(ps_fam))
colnames(asv_table) <- tax_table_df$Family
```

---

## 2. Modelagem com MaAsLin3

Agora que temos nossa fatia de dados contendo apenas mel e líquen, podemos rodar o modelo linear para perguntar: *Quais famílias microbianas estão significativamente enriquecidas ou depletadas no mel em relação ao líquen?*

Adicione o restante do código ao seu script:

```r
# 6. Rodar a modelagem linear
cat("\nExecutando MaAsLin3...\n")
fit_data <- maaslin3(
    input_data = asv_table,
    input_metadata = metadata,
    output = "results/maaslin3_sliced",  # Diretório de saída exclusivo
    normalization = "TSS",               # Total Sum Scaling (fração)
    transform = "LOG",                   # Transformação Log
    fixed_effects = "group",             # Queremos ver o efeito do 'group'
    reference = "group,liquen"           # O líquen será a nossa linha de base!
)

cat("\nAnálise MaAsLin3 concluída! Resultados em results/maaslin3_sliced/\n")
```

> [!TIP]
> **Por que `reference = "group,liquen"`?**
> A referência define o "nível zero" do modelo. Ao setar o líquen como referência, os coeficientes estatísticos que o programa calcular para a categoria **mel** nos dirão se a abundância no mel aumentou (coeficiente positivo) ou diminuiu (coeficiente negativo) em relação ao que existe no líquen.

---

## 3. Submissão ao SLURM

Como este é um dataset pequeno ("slice"), a análise será muito rápida, mas é uma boa prática rodarmos tudo pelo cluster. Crie o arquivo `scripts/03_maaslin_sliced.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=maaslin_slice
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/maaslin_slice_%j.out
#SBATCH --error=logs/maaslin_slice_%j.err

mkdir -p logs
mkdir -p results/maaslin3_sliced

# Carregar o R e executar o script
module load HPC/R/4.3.3
Rscript scripts/maaslin3_sliced.R
```

Submeta com:
```bash
sbatch scripts/03_maaslin_sliced.sh
```

---

## 4. Interpretando os Resultados do Slice

Quando o job terminar, navegue até a pasta `results/maaslin3_sliced`. Você encontrará:

### 📄 `all_results.tsv`
Abra esta tabela. Observe as colunas:
- `feature`: Nome da Família bacteriana.
- `coef`: Coeficiente de abundância. Se for positivo (ex: `1.5`), significa que a abundância desta família é **maior no mel** do que no líquen. Se for negativo, é **menor**.
- `qval`: p-valor ajustado. Se for `< 0.05`, a diferença é estatisticamente significativa!

### 📊 Pasta `figures/`
O programa gerará automaticamente:
- **Heatmap (`heatmap.pdf`)**: Um gráfico de calor mostrando as famílias mais significativas e a direção da associação com a amostra de "mel".
- **Scatterplots / Boxplots**: PDFs individuais para cada família que deu diferença estatística significativa, mostrando claramente a abundância nos dois grupos.

---
**Próximo Passo:** Após entender a mecânica do modelo neste cenário simples, você está pronto para integrar múltiplos bancos de dados e adicionar correção de covariáveis (efeito de lote) com o [Tutorial Principal da Parte 03](TUTORIAL.md).
