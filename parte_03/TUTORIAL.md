# Tutorial: Integração e Análise de Abundância Diferencial com MaAsLin3

> **Projeto:** LZT0693 — Iniciação Científica em Biotecnologia  
> **Objetivo:** Combinar dados brutos das Partes 01 e 02, corrigir efeitos de lote (batch effects) e executar análise de associação estatística multivariada.

---

## Sumário

1. [Contexto: Por que integrar dados de diferentes aulas?](#1-contexto-por-que-integrar-dados-de-diferentes-aulas)
2. [Visão geral do fluxo de análise](#2-visão-geral-do-fluxo-de-análise)
3. [Pré-requisitos](#3-pré-requisitos)
4. [Passo a passo: Construindo o script R](#4-passo-a-passo-construindo-o-script-r)
   - 4.1 [Configurar ambiente e caminhos](#41-configurar-ambiente-e-caminhos)
   - 4.2 [Carregar tabelas de ASVs e taxonomia](#42-carregar-tabelas-de-asvs-e-taxonomia)
   - 4.3 [Mesclar tabelas de ASVs](#43-mesclar-tabelas-de-asvs)
   - 4.4 [Combinar taxonomias](#44-combinar-taxonomias)
   - 4.5 [Preparar e limpar metadados](#45-preparar-e-limpar-metadados)
   - 4.6 [Calcular profundidade de sequenciamento (`reads`)](#46-calcular-profundidade-de-sequenciamento-reads)
   - 4.7 [Agrupar por nível taxonômico (`tax_glom`)](#47-agrupar-por-nível-taxonômico-tax_glom)
   - 4.8 [Executar o MaAsLin3](#48-executar-o-maaslin3)
5. [Interpretando os resultados](#5-interpretando-os-resultados)
6. [Executando o script completo](#6-executando-o-script-completo)
7. [Resolução de problemas (Troubleshooting)](#7-resolução-de-problemas-troubleshooting)
8. [Glossário de termos](#8-glossário-de-termos)
9. [Referências bibliográficas](#9-referências-bibliográficas)

---

## 1. Contexto: Por que integrar dados de diferentes aulas?

Nas aulas anteriores, vocês processaram dois conjuntos de amostras de forma independente:

| | Parte 01 | Parte 02 |
|---|---|---|
| **Amostras** | Líquen, musgo, mel | Grupo A, Grupo B |
| **Controle negativo** | `S813_16Sv3v4_NN` | `S813_NN_V3V4` |
| **Sequenciamento** | Corrida/lote 1 | Corrida/lote 2 |

Cada aula gerou sua própria tabela de ASVs (`seqtab_nochim.rds`) e classificação taxonômica (`taxa.rds`). Agora, queremos **combinar tudo** em uma análise única.

> **⚠️ O problema:** Como as amostras foram sequenciadas em momentos diferentes, existe uma variação técnica chamada **efeito de lote** (*batch effect*). Se não corrigirmos, o modelo pode confundir diferenças biológicas reais com artefatos técnicos do sequenciamento.

> **✅ A solução:** Usar o **MaAsLin3** (*Microbiome Multivariable Association with Linear Models*), que permite modelar a abundância de cada táxon incluindo covariáveis de ajuste como `batch` (lote) e `reads` (profundidade de sequenciamento). Assim, a fórmula estatística fica:
>
> $$\text{Abundância} \sim \text{group} + \text{batch} + \text{reads}$$

---

## 2. Visão geral do fluxo de análise

Antes de mergulhar no código, veja o panorama do que faremos:

```
┌─────────────────┐    ┌─────────────────┐
│    Parte 01     │    │    Parte 02     │
│ seqtab_nochim   │    │ seqtab_nochim   │
│ taxa.rds        │    │ taxa.rds        │
│ metadata.csv    │    │ metadata.csv    │
└────────┬────────┘    └────────┬────────┘
         │                      │
         └──────────┬───────────┘
                    ▼
        ┌───────────────────────┐
        │  Mesclar ASVs (DADA2) │  ← mergeSequenceTables()
        └───────────┬───────────┘
                    ▼
        ┌───────────────────────┐
        │  Combinar taxonomias  │  ← rbind + deduplicate
        └───────────┬───────────┘
                    ▼
        ┌───────────────────────┐
        │  Harmonizar metadados │  ← adicionar coluna "batch"
        │  + remover controles  │     + coluna "reads"
        └───────────┬───────────┘
                    ▼
        ┌───────────────────────┐
        │  Agrupar por nível    │  ← tax_glom("Family")
        │  taxonômico (Phyloseq)│
        └───────────┬───────────┘
                    ▼
        ┌───────────────────────┐
        │       MaAsLin3        │  ← ~ group + batch + reads
        │  Abundância diferencial│
        └───────────┬───────────┘
                    ▼
           results/maaslin3/
```

---

## 3. Pré-requisitos

Antes de começar, verifique se:

- [x] Você concluiu os pipelines DADA2 da **Parte 01** e da **Parte 02**.
- [x] Os arquivos `seqtab_nochim.rds` e `taxa.rds` existem nos diretórios `results/dada2/` de ambas as partes.
- [x] Os arquivos `metadata.csv` estão corretos em ambas as partes.
- [x] O pacote `maaslin3` está instalado no R (se não, rode: `BiocManager::install("maaslin3")`).

---

## 4. Passo a passo: Construindo o script R

Crie o arquivo `scripts/maaslin_analysis.R` dentro da pasta `parte_03/` e cole cada bloco de código a seguir. Você pode colar tudo de uma vez ou ir colando passo a passo no console do RStudio para acompanhar cada resultado.

### 4.1 Configurar ambiente e caminhos

Carregamos os pacotes e definimos os caminhos para os dados das etapas anteriores.

> **📝 Importante:** Substitua os caminhos genéricos (`caminho/para/...`) pelos caminhos reais no seu sistema. Por exemplo, se você clonou o repositório em `/home/seu_usuario/LZT0693/`, os caminhos seriam algo como:
> ```r
> seqtab1_path <- "/home/seu_usuario/LZT0693/parte_01/results/dada2/seqtab_nochim.rds"
> ```

```r
# ==============================================================================
# Análise de Abundância Diferencial com MaAsLin3 (Integração Parte 01 & 02)
# Disciplina LZT0693 - Iniciação Científica em Biotecnologia
# ==============================================================================

# Carregar bibliotecas necessárias
library(dada2)
library(maaslin3)
library(dplyr)
library(phyloseq)

# ── Caminhos das tabelas de ASVs ──
seqtab1_path <- "caminho/para/parte_01/seqtab_nochim.rds"
seqtab2_path <- "caminho/para/parte_02/seqtab_nochim.rds"

# ── Caminhos dos metadados ──
meta1_path <- "caminho/para/parte_01/metadata.csv"
meta2_path <- "caminho/para/parte_02/metadata.csv"

# ── Caminhos das classificações taxonômicas ──
taxa1_path <- "caminho/para/parte_01/taxa.rds"
taxa2_path <- "caminho/para/parte_02/taxa.rds"

# ── Diretório de saída ──
output_dir <- "../results/maaslin3"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
```

### 4.2 Carregar tabelas de ASVs e taxonomia

Fazemos a leitura dos arquivos `.rds` salvos pelo DADA2 nas etapas anteriores. O script verifica se os arquivos existem antes de prosseguir.

```r
# Verificar se os arquivos existem
if (!file.exists(seqtab1_path) || !file.exists(seqtab2_path) ||
    !file.exists(taxa1_path) || !file.exists(taxa2_path)) {
  stop("⚠️ Erro: Arquivos não encontrados. Verifique os caminhos definidos no Passo 4.1!")
}

# Carregar tabelas de ASVs (matriz: amostras × sequências)
seqtab1 <- readRDS(seqtab1_path)
seqtab2 <- readRDS(seqtab2_path)

# Carregar classificações taxonômicas (matriz: sequências × níveis taxonômicos)
tax1 <- readRDS(taxa1_path)
tax2 <- readRDS(taxa2_path)

cat("Parte 01:", nrow(seqtab1), "amostras,", ncol(seqtab1), "ASVs\n")
cat("Parte 02:", nrow(seqtab2), "amostras,", ncol(seqtab2), "ASVs\n")
```

> **🔍 O que esperar:** Você deverá ver algo como:
> ```
> Parte 01: 8 amostras, 1523 ASVs
> Parte 02: 17 amostras, 2841 ASVs
> ```
> Os números exatos dependerão dos seus dados.

### 4.3 Mesclar tabelas de ASVs

A função `mergeSequenceTables()` do DADA2 compara as sequências de DNA de cada ASV e **une automaticamente** aquelas que são idênticas entre os dois lotes. ASVs exclusivas de um lote são mantidas com contagem zero nas amostras do outro.

```r
# Mesclar as tabelas (ASVs com sequências idênticas são pareadas)
merged_seqtab <- mergeSequenceTables(seqtab1, seqtab2)
cat("Tabela mesclada:", nrow(merged_seqtab), "amostras,", ncol(merged_seqtab), "ASVs\n")
```

> **🔍 O que esperar:** O número de ASVs na tabela mesclada será **menor** que a soma das duas tabelas individuais, pois ASVs idênticas são unificadas.

### 4.4 Combinar taxonomias

Unimos as tabelas de classificação taxonômica e removemos as linhas duplicadas (ASVs que apareciam em ambos os lotes terão a mesma sequência como nome de linha).

```r
# Combinar taxonomias e remover duplicatas
combined_tax <- rbind(tax1, tax2)
combined_tax <- combined_tax[!duplicated(rownames(combined_tax)), ]
```

### 4.5 Preparar e limpar metadados

Aqui fazemos três coisas importantes:
1. **Adicionamos a coluna `batch`** — identifica de qual aula/corrida cada amostra veio.
2. **Removemos os controles negativos** — eles não devem entrar na modelagem estatística.
3. **Sincronizamos** as tabelas — garantimos que as mesmas amostras estejam na tabela de ASVs e nos metadados.

```r
# Carregar metadados
meta1 <- read.csv(meta1_path, row.names = 1)
meta2 <- read.csv(meta2_path, row.names = 1)

# Adicionar coluna indicando o lote de origem
meta1$batch <- "parte_01"
meta2$batch <- "parte_02"

# Unir os metadados em uma única tabela
merged_metadata <- rbind(meta1, meta2)

# Remover controles negativos
control_samples <- c("S813_16Sv3v4_NN", "S813_NN_V3V4")
merged_metadata <- merged_metadata[!rownames(merged_metadata) %in% control_samples, ]

# Manter apenas amostras presentes em ambas as tabelas
samples_to_keep <- intersect(rownames(merged_metadata), rownames(merged_seqtab))
merged_seqtab <- merged_seqtab[samples_to_keep, ]
merged_metadata <- merged_metadata[samples_to_keep, ]
```

### 4.6 Calcular profundidade de sequenciamento (`reads`)

A profundidade de sequenciamento (total de reads por amostra) pode variar bastante entre amostras e entre lotes. Essa variação é uma das principais fontes de confundimento em estudos de microbioma: uma amostra que foi sequenciada com maior profundidade de biblioteca terá **mais reads mapeados** e, consequentemente, mais chance de detectar ASVs raras — mesmo que a comunidade microbiana subjacente seja idêntica.

Isso significa que diferenças aparentes no número de ASVs entre duas amostras podem não refletir diferenças biológicas reais, mas simplesmente o fato de que uma amostra recebeu mais sequências do que a outra durante o sequenciamento. Se não corrigirmos, o modelo pode interpretar esses artefatos de cobertura como sinais biológicos.

A inclusão da profundidade de sequenciamento como covariável é uma prática recomendada pelos próprios desenvolvedores do MaAsLin3 (equipe do Huttenhower Lab, Harvard T.H. Chan School of Public Health). No [tutorial oficial do MaAsLin3](https://github.com/biobakery/maaslin3) e no artigo de referência ([Mallick et al., 2024 — *Nature Methods*](https://doi.org/10.1038/s41592-024-02394-8)), os autores orientam incluir variáveis técnicas conhecidas (como lote e profundidade) como covariáveis no modelo para isolar corretamente os efeitos biológicos de interesse.

Adicionamos essa informação aos metadados como uma coluna numérica:

```r
# Calcular o total de reads por amostra e adicionar aos metadados
merged_metadata$reads <- rowSums(merged_seqtab)
```

> **💡 Resumindo:** Incluir `reads` na fórmula garante que o MaAsLin3 desconte as diferenças de profundidade de sequenciamento. Sem essa covariável, uma amostra com 50.000 reads parecerá "mais diversa" do que uma com 5.000 reads, e o modelo pode confundir esse artefato com uma diferença biológica real entre os grupos.

### 4.6b Configurar variáveis categóricas (Fatores) e Nível de Referência (Baselines)

O R e o MaAsLin3 precisam entender quais colunas são categóricas (fatores) e qual é o grupo "controle" ou linha de base (baseline) para a comparação. 

1. **Fatores:** Se uma variável contendo texto (ex: `"rupicula"`, `"epifita"`) não for explicitamente convertida em fator, o R pode lê-la como `character`. Algumas bibliotecas a interpretam de forma contínua, gerando gráficos de tendências lineares estranhos em vez de boxplots por grupo.
2. **Nível de Referência:** Por padrão, o R ordena os níveis de forma alfabética (no caso do tratamento: `epifita` seria a referência). Para alterar o controle para `rupicula` (ou seja, avaliar `epifita vs rupicula` e `terricola vs rupicula`), usamos `relevel()`.

```r
# Converter variáveis categóricas para fator
merged_metadata$treatment <- factor(merged_metadata$treatment)
merged_metadata$group     <- factor(merged_metadata$group)
merged_metadata$batch     <- factor(merged_metadata$batch)

# Definir "rupicula" como nível de referência para a comparação de tratamentos
if ("rupicula" %in% levels(merged_metadata$treatment)) {
  merged_metadata$treatment <- relevel(merged_metadata$treatment, ref = "rupicula")
}
```


### 4.7 Agrupar por nível taxonômico (`tax_glom`)

Até aqui, cada coluna da nossa tabela representa uma ASV individual (uma sequência de DNA única). São centenas ou milhares de colunas, muitas delas com contagens muito baixas (*esparsas*).

Para facilitar a interpretação biológica e aumentar o poder estatístico, podemos **agrupar** as ASVs que pertencem à mesma família (ou outro nível taxonômico), somando suas contagens.

> **📋 Níveis taxonômicos disponíveis para agrupamento:**
>
> | Valor em R | Nível | Descrição |
> |---|---|---|
> | `"Kingdom"` | Reino | Muito amplo (Bacteria, Archaea) |
> | `"Phylum"` | Filo | Agrupamento amplo (ex: Firmicutes, Proteobacteria) |
> | `"Class"` | Classe | Intermediário |
> | `"Order"` | Ordem | Intermediário |
> | **`"Family"`** | **Família** | **← Usado neste tutorial (bom equilíbrio entre resolução e robustez)** |
> | `"Genus"` | Gênero | Alta resolução, mas pode ser mais esparso |
> | `"Species"` | Espécie | Muito específico (muitos NAs na classificação) |

```r
# Criar objeto Phyloseq
ps <- phyloseq(
  otu_table(merged_seqtab, taxa_are_rows = FALSE),
  sample_data(merged_metadata),
  tax_table(combined_tax)
)

# Agrupar as ASVs no nível de Família
# Para outro nível, altere "Family" para o desejado (ver tabela acima)
ps_glom <- tax_glom(ps, taxrank = "Family")

# Renomear os táxons agrupados
# Por padrão, o tax_glom mantém a sequência de DNA representante da ASV como nome de coluna.
# Aqui extraímos a taxonomia correspondente e criamos nomes amigáveis no formato 'Phylum__Family'.
tax_mat <- as(tax_table(ps_glom), "matrix")
new_names <- ifelse(
  !is.na(tax_mat[, "Family"]),
  paste(tax_mat[, "Phylum"], tax_mat[, "Family"], sep = "__"),
  ifelse(
    !is.na(tax_mat[, "Order"]),
    paste(tax_mat[, "Phylum"], tax_mat[, "Order"], "unclassified", sep = "__"),
    paste(tax_mat[, "Phylum"], "unclassified", sep = "__")
  )
)
# Garante que os nomes sejam únicos e atribui de volta ao objeto phyloseq
taxa_names(ps_glom) <- make.unique(new_names, sep = "_")

# Extrair os dados finais limpos para usar no MaAsLin3
final_seqtab   <- as(otu_table(ps_glom), "matrix")
final_metadata <- as(sample_data(ps_glom), "data.frame")

cat("Amostras para análise:", nrow(final_metadata), "\n")
cat("Famílias taxonômicas:", ncol(final_seqtab), "\n")
```

> **🔍 O que esperar:** O número de colunas (features) cairá drasticamente. Por exemplo, de ~3.000 ASVs para ~150 famílias. Agora, os nomes dos táxons serão legíveis (ex: `Firmicutes__Lactobacillaceae`), facilitando a análise posterior.

### 4.8 Executar o MaAsLin3 com seleção dinâmica de covariáveis

Finalmente, rodamos o modelo. A fórmula básica desejada seria `~ group + treatment + batch + reads`. Porém, se filtrarmos os dados (ex: selecionando apenas um grupo), algumas dessas variáveis passam a ter **apenas 1 valor/nível**, o que causará o erro clássico de contrastes no R. 

Para tornar o script robusto, adicionamos uma etapa que remove automaticamente da fórmula qualquer preditor sem variância real:

```r
# Definir preditores candidatos
formula_vars <- c("group", "treatment", "batch", "reads")

# Filtrar para manter apenas variáveis com variabilidade (fatores com >= 2 níveis ou numéricos)
vars_to_keep <- sapply(formula_vars, function(v) {
  if (!v %in% colnames(final_metadata)) return(FALSE)
  col <- final_metadata[[v]]
  if (is.factor(col) || is.character(col)) {
    return(length(unique(na.omit(col))) >= 2)
  }
  return(TRUE)
})

vars_kept    <- formula_vars[vars_to_keep]
vars_dropped <- formula_vars[!vars_to_keep]

if (length(vars_dropped) > 0) {
  cat("⚠️ Variáveis removidas da fórmula (apenas 1 nível):", paste(vars_dropped, collapse = ", "), "\n")
}

# Criar a string de fórmula final
maaslin_formula <- paste("~", paste(vars_kept, collapse = " + "))
cat("✅ Fórmula final construída:", maaslin_formula, "\n")

# Executar o MaAsLin3
maaslin_results <- maaslin3(
  input_data       = final_seqtab,
  input_metadata   = final_metadata,
  output           = output_dir,
  formula          = maaslin_formula,
  normalization    = "TSS",
  transform        = "LOG",
  min_abundance    = 0.0001,
  min_prevalence   = 0.1,
  max_significance = 0.05
)

cat("\n>> Os resultados foram salvos no diretório:", output_dir, "\n")
```


#### Detalhamento dos parâmetros

| Parâmetro | Valor usado | O que faz | Alternativas |
|-----------|-------------|-----------|--------------|
| `input_data` | `merged_seqtab_glom` | Matriz de abundâncias (amostras × táxons). Aqui usamos os dados já agrupados por família. | Pode ser a tabela de ASVs original (`merged_seqtab`) se não quiser agrupar. |
| `input_metadata` | `final_metadata` | Tabela com as variáveis clínicas/experimentais de cada amostra. | — |
| `output` | `output_dir` | Pasta onde os resultados serão salvos (tabelas `.tsv`, gráficos `.pdf`). | Qualquer caminho válido. |
| `formula` | `"~ group + batch + reads"` | Define o modelo estatístico. `group` é a variável de interesse; `batch` e `reads` são covariáveis de ajuste. | Remover `reads` se a profundidade for homogênea; adicionar outras variáveis como `treatment` ou `sex`. |
| `normalization` | `"TSS"` | **Total Sum Scaling**: divide cada contagem pelo total da amostra, convertendo em proporções (abundância relativa). É o método mais intuitivo e recomendado para dados composicionais de microbioma. | `"CLR"` (Centered Log-Ratio, mais robusto para composicionalidade mas menos intuitivo), `"CSS"` (Cumulative Sum Scaling), `"NONE"` (sem normalização). |
| `transform` | `"LOG"` | Aplica transformação logarítmica (`log(x + pseudocount)`) após a normalização. Estabiliza a variância e reduz o peso de táxons hiper-abundantes, tornando a distribuição mais adequada para modelos lineares. | `"AST"` (Arcsine Square Root, alternativa para dados proporcionais), `"NONE"` (sem transformação). |
| `min_abundance` | `0.0001` | Abundância relativa mínima para que um táxon seja considerado na análise. Táxons com abundância média abaixo deste limiar são descartados antes da modelagem. | Aumentar para `0.001` se quiser focar apenas em táxons mais abundantes; diminuir para `0.00001` se quiser incluir mais táxons raros. |
| `min_prevalence` | `0.1` | Prevalência mínima: o táxon deve estar presente em pelo menos 10% das amostras. Remove táxons muito esporádicos que geram ruído estatístico. | Aumentar para `0.2` ou `0.3` para análises mais conservadoras; diminuir para `0.05` se tiver poucas amostras. |
| `max_significance` | `0.05` | Limiar de significância para o q-value (p-valor ajustado por FDR). Apenas resultados com q-value abaixo deste valor serão reportados como significativos. | `0.1` para análises exploratórias; `0.01` para critérios mais rigorosos. |
| `plot_heatmap` | `TRUE` | Gera um heatmap (`.pdf`) com os coeficientes de associação dos táxons significativos. | `FALSE` para desativar. |
| `plot_scatter` | `TRUE` | Gera gráficos de dispersão/boxplot individuais para cada táxon significativo. | `FALSE` para desativar (útil se houver muitos resultados significativos). |

---

## 5. Interpretando os resultados

Os resultados serão gravados em `results/maaslin3/` e contêm:

### 📄 Tabela principal: `all_results.tsv`

| Coluna | O que significa |
|--------|-----------------|
| `feature` | Nome do táxon (família, gênero, etc.) |
| `metadata` | A covariável testada (`group`, `batch` ou `reads`) |
| `coef` | Coeficiente do modelo: **positivo** = mais abundante no grupo testado; **negativo** = menos abundante |
| `pval` | p-valor bruto |
| `qval` | p-valor corrigido (FDR). **Táxons com `qval < 0.05` são considerados estatisticamente significativos** |

### 📊 Visualizações

- **Heatmap (`maaslin3_heatmap.pdf`)**: Mostra os coeficientes de associação para todos os táxons significativos de forma consolidada.
- **Scatter plots / Boxplots**: Gráficos individuais para cada táxon significativo, mostrando a distribuição de abundância relativa entre os grupos.

> **🔍 O que procurar:** Foque nas linhas onde `metadata == "group"` e `qval < 0.05`. Essas são as famílias bacterianas que apresentam diferença estatisticamente significativa entre seus grupos biológicos, **já corrigidas** para efeito de lote e profundidade de sequenciamento.

### 🧪 Exemplo concreto de interpretação

Suponha que você encontre a seguinte linha no arquivo `all_results.tsv`:

| feature | metadata | value | coef | stderr | pval | qval |
|---------|----------|-------|------|--------|------|------|
| Lachnospiraceae | group | musgo | 1.82 | 0.41 | 0.0003 | 0.012 |

Como interpretar:

1. **`feature = Lachnospiraceae`**: A família bacteriana analisada.
2. **`metadata = group`**: O resultado se refere à variável de interesse (grupo biológico), não ao lote ou profundidade.
3. **`value = musgo`**: O grupo específico sendo comparado (em relação ao grupo de referência, que o R define alfabeticamente — neste caso, seria `liquen`).
4. **`coef = 1.82`**: O coeficiente é **positivo**, indicando que a família Lachnospiraceae é **mais abundante** nas amostras de musgo do que nas de líquen (grupo de referência). Em escala logarítmica, isso corresponde a aproximadamente `10^1.82 ≈ 66 vezes` mais abundante.
5. **`qval = 0.012`**: O q-value (p-valor corrigido para múltiplos testes via FDR) é **menor que 0.05**, portanto este resultado é **estatisticamente significativo**.

> **💡 Atenção ao grupo de referência:** O R ordena os níveis de fatores alfabeticamente por padrão. No exemplo acima com grupos `liquen`, `mel` e `musgo`, o grupo de referência seria `liquen`. Os coeficientes de `mel` e `musgo` são sempre relativos a ele. Para alterar o grupo de referência, use `relevel()` antes de rodar o MaAsLin3:
> ```r
> final_metadata$group <- relevel(factor(final_metadata$group), ref = "mel")
> ```

---

## 6. Executando o script completo

Após construir o script (ou colar todos os blocos acima no arquivo `scripts/maaslin_analysis.R`):

1. **Edite os caminhos** no Passo 4.1 para apontar para os seus arquivos reais.
2. **Verifique os metadados** de ambas as partes (`metadata.csv`).
3. **Execute** de uma das formas:
   - No **console do RStudio**: copie e cole cada bloco sequencialmente.
   - No **terminal**: rode o script inteiro com:
     ```bash
     Rscript scripts/maaslin_analysis.R
     ```

---

## 7. Resolução de problemas (Troubleshooting)

Abaixo estão os erros mais comuns e como resolvê-los:

### ❌ `Error in readRDS(...): cannot open connection`

**Causa:** O caminho para o arquivo `.rds` está incorreto ou o arquivo não existe.

**Solução:**
- Verifique se você editou corretamente os caminhos no Passo 4.1.
- Confirme que os pipelines das Partes 01 e 02 foram executados com sucesso.
- Use `file.exists("seu/caminho/aqui.rds")` no console do R para testar se o arquivo é encontrado.

### ❌ `Error in rbind(...): number of columns of arguments do not match`

**Causa:** As tabelas de metadados (`metadata.csv`) das duas partes têm colunas diferentes.

**Solução:**
- Abra ambos os arquivos CSV e verifique se as colunas são idênticas (mesmo nome e mesma ordem). As colunas esperadas são: `sample_name`, `group`, `treatment`, `replicate`, `description`.
- Se uma das partes tiver colunas extras, remova-as ou adicione-as à outra parte com valores `NA`.

### ❌ `Error in tax_glom(...): taxrank must be a column name of the taxonomy table`

**Causa:** O nível taxonômico especificado (ex: `"Family"`) não existe na tabela de taxonomia.

**Solução:**
- Verifique os nomes das colunas da sua tabela taxonômica com:
  ```r
  colnames(tax_table(ps))
  ```
- Os nomes devem ser: `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species`. Se estiverem diferentes (ex: `family` em minúsculo), ajuste o parâmetro `taxrank` de acordo.

### ❌ `Error in maaslin3(...): could not find function "maaslin3"`

**Causa:** O pacote `maaslin3` não está instalado.

**Solução:**
```r
BiocManager::install("maaslin3")
library(maaslin3)
```

### ❌ Nenhum resultado significativo (`all_results.tsv` vazio ou sem `qval < 0.05`)

**Causa:** Pode ser falta de poder estatístico (poucas amostras), filtros muito restritivos, ou simplesmente não há diferença biológica detectável.

**O que tentar:**
- Reduza `min_prevalence` para `0.05` para incluir mais táxons.
- Reduza `min_abundance` para `0.00001`.
- Aumente `max_significance` para `0.1` (análise exploratória).
- Tente agrupar em um nível taxonômico mais alto (ex: `"Order"` ou `"Phylum"` ao invés de `"Family"`).

### ❌ `Error in contrasts<-` (`contrasts can be applied only to factors with 2 or more levels`)

**Causa:** Uma variável categórica na fórmula tem apenas **1 nível** nos dados após os filtros aplicados (por exemplo, você filtrou apenas pelo grupo `"musgo"`, deixando a variável `group` com apenas um nível, ou todas as amostras são da `"parte_02"`, deixando `batch` com apenas um nível). O modelo linear de regressão exige variação (pelo menos dois grupos/valores diferentes) para calcular contrastes e estimar coeficientes.

**Solução:**
- A lógica de seleção de covariáveis adicionada no script na etapa 8a já faz esse descarte dinâmico automaticamente, excluindo preditores constantes e mostrando alertas como: `⚠️ Variáveis removidas da fórmula (apenas 1 nível): group`.
- Se você quiser realizar a comparação entre tratamentos para o grupo `"musgo"`, a fórmula gerada automaticamente será `~ treatment + batch + reads`. Não é necessário alterar a fórmula manualmente!

### ⚠️ O gráfico gerou tendências lineares (linhas) para variáveis de grupos discretos

**Causa:** A coluna de agrupamento experimental (ex: `treatment`) foi lida como `character` ou `numeric`, fazendo com que o MaAsLin3 assumisse que a variável é contínua e realizasse uma regressão com linha de tendência.

**Solução:**
- Certifique-se de que a etapa **6b** de conversão para fatores está habilitada no script:
  ```r
  merged_metadata$treatment <- factor(merged_metadata$treatment)
  ```
  Isso instrui o R a tratar os valores como grupos categóricos, forçando o MaAsLin3 a desenhar boxplots/beeswarms em vez de linhas.

### ⚠️ O script roda mas os resultados parecem estranhos

**Verificações recomendadas:**
- Confirme que os nomes das amostras nos metadados batem exatamente com os nomes das linhas na tabela de ASVs (`rownames`). Diferenças de maiúsculas/minúsculas ou espaços extras causam problemas silenciosos.
- Verifique se a coluna `group` nos metadados contém os valores corretos (sem typos).
- Inspecione a tabela de metadados mesclada com `head(merged_metadata)` e `str(merged_metadata)` para confirmar que `batch` e `reads` foram adicionados corretamente.


---

## 8. Glossário de termos

| Termo | Significado |
|-------|-------------|
| **ASV** | *Amplicon Sequence Variant*. Sequência de DNA única inferida pelo DADA2 após remoção de erros de sequenciamento. Diferente de OTUs, as ASVs têm resolução de nucleotídeo único. |
| **Batch effect** | Variação técnica sistemática introduzida por diferenças entre corridas de sequenciamento, preparação de biblioteca ou outros fatores não biológicos. |
| **Covariável** | Variável incluída no modelo estatístico para "controlar" seu efeito. No nosso caso, `batch` e `reads` são covariáveis — não queremos testá-las, mas sim descontar sua influência. |
| **DADA2** | Pipeline de bioinformática que infere ASVs a partir de dados de amplicon, corrigindo erros de sequenciamento com um modelo probabilístico. |
| **FDR** | *False Discovery Rate*. Método de correção para múltiplos testes que controla a proporção esperada de falsos positivos entre os resultados significativos. O q-value é o p-valor ajustado por FDR. |
| **Library size** | Número total de reads (sequências) obtido para uma amostra. Sinônimo de "profundidade de sequenciamento". |
| **MaAsLin3** | *Microbiome Multivariable Association with Linear Models*. Ferramenta que modela a abundância de cada táxon contra variáveis de interesse usando modelos lineares generalizados. |
| **Normalização** | Processo de tornar as contagens comparáveis entre amostras. TSS (Total Sum Scaling) divide cada contagem pelo total da amostra, gerando proporções. |
| **Phyloseq** | Pacote do R para manipulação e análise de dados de microbioma. Permite agrupar ASVs, calcular diversidade, e gerar visualizações. |
| **p-valor** | Probabilidade de observar um resultado tão extremo quanto o obtido, assumindo que a hipótese nula ("não há diferença") é verdadeira. Quanto menor, mais evidência contra a hipótese nula. |
| **q-valor** | p-valor corrigido para múltiplos testes (via FDR). Quando testamos centenas de táxons simultaneamente, o q-value evita que encontremos diferenças "significativas" por puro acaso. |
| **tax_glom** | *Taxonomic Glomming*. Função do Phyloseq que agrega ASVs pelo nível taxonômico desejado (ex: Família), somando as contagens de todas as ASVs classificadas naquele grupo. |
| **TSS** | *Total Sum Scaling*. Método de normalização que divide cada contagem pelo total de reads da amostra, convertendo para abundância relativa (proporções que somam 1). |
| **Transformação LOG** | Aplicação de logaritmo nas abundâncias normalizadas. Reduz a influência de táxons hiper-abundantes e torna a distribuição dos dados mais simétrica, adequada para modelos lineares. |

---

## 9. Referências bibliográficas

1. **MaAsLin3** — Mallick, H. et al. (2024). Multivariable association discovery in population-scale meta-omics studies. *Nature Methods*. DOI: [10.1038/s41592-024-02394-8](https://doi.org/10.1038/s41592-024-02394-8)
2. **DADA2** — Callahan, B.J. et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583. DOI: [10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
3. **Phyloseq** — McMurdie, P.J. & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE*, 8(4), e61217. DOI: [10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)
4. **SILVA** — Quast, C. et al. (2013). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Research*, 41, D590–D596. DOI: [10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
5. **Tutorial oficial MaAsLin3** — Huttenhower Lab, Harvard T.H. Chan School of Public Health. Disponível em: [https://github.com/biobakery/maaslin3](https://github.com/biobakery/maaslin3)
