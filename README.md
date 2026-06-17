# LZT0693 — Análise de Amplicon 16S rRNA (V3-V4)

Este repositório reúne os tutoriais, scripts e materiais didáticos para a análise de sequenciamento de amplicons da região **16S rRNA (V3-V4)** desenvolvido na disciplina **LZT0693 — Iniciação Científica em Biotecnologia**.

O projeto está dividido em três módulos complementares, cobrindo desde o controle de qualidade de dados brutos até a modelagem estatística multivariada com correção de efeito de lote.

---

## 📂 Estrutura Geral do Repositório

```
LZT0693/
├── README.md               # ← Este arquivo (visão geral)
├── parte_01/               # Módulo 1 — Pipeline DADA2 introdutório
│   ├── TUTORIAL.md
│   ├── metadata.csv
│   ├── scripts/
│   └── slides/
├── parte_02/               # Módulo 2 — Pipeline com novos tratamentos
│   ├── TUTORIAL.md
│   ├── metadata.csv
│   └── scripts/
└── parte_03/               # Módulo 3 — Integração e MaAsLin3
    ├── TUTORIAL.md
    └── scripts/
```

---

## 🔬 Descrição dos Módulos

### Parte 01 — Introdução ao DADA2 e Ecologia Microbiana

Foco no entendimento do fluxo de qualidade de sequências brutas e geração de tabelas taxonômicas.

* **Amostras**: 8 bibliotecas (líquen, musgo, mel e 1 controle negativo).
* **Etapas**: FastQC → Cutadapt → DADA2 (filtragem, denoising, quimeras, taxonomia SILVA) → Phyloseq (diversidade alfa/beta).
* **Tutorial completo**: [`parte_01/TUTORIAL.md`](parte_01/TUTORIAL.md)

### Parte 02 — Pipeline de Produção e Novos Grupos

Aplicação do mesmo pipeline em um conjunto de amostras maior e com desenho experimental de dois grupos.

* **Amostras**: 17 bibliotecas (12 do Grupo A, 4 do Grupo B e 1 controle negativo).
* **Etapas**: Idênticas à Parte 01, com foco na comparação `Grupo_A` vs. `Grupo_B`.
* **Tutorial completo**: [`parte_02/TUTORIAL.md`](parte_02/TUTORIAL.md)

### Parte 03 — Integração e Modelagem Estatística com MaAsLin3

Integração de múltiplos conjuntos de dados e análise de abundância diferencial multivariada.

* **Objetivo**: Combinar as tabelas de ASVs e taxonomias das Partes 01 e 02, agrupar por nível taxonômico (ex. Família) e rodar modelagens lineares com correção de efeito de lote e profundidade de sequenciamento.
* **Etapas**: Mesclagem de ASVs → Harmonização de metadados → `tax_glom` (Phyloseq) → MaAsLin3.
* **Tutorial completo**: [`parte_03/TUTORIAL.md`](parte_03/TUTORIAL.md)

---

## ⚙️ Configuração Específica por Módulo

Cada módulo precisa de configurações específicas **antes de iniciar** a análise. A tabela abaixo resume os passos de configuração obrigatórios:

### Parte 01 e Parte 02 (Pipeline DADA2)

Ambas as partes seguem os mesmos pré-requisitos. Execute os comandos abaixo de **dentro do diretório** de cada parte (`cd parte_01` ou `cd parte_02`):

**1. Criar os links simbólicos para dados e bancos de dados:**

```bash
# Parte 01
cd parte_01
ln -s /home/cursos/LCoutinho202604/data_shared/aula_02/data     data
ln -s /home/cursos/LCoutinho202604/data_shared/aula_02/databases databases

# Parte 02
cd ../parte_02
ln -s /home/cursos/LCoutinho202604/data_shared/aula_03/data     data
ln -s /home/cursos/LCoutinho202604/data_shared/aula_03/databases databases
```

**2. Carregar os módulos de software no terminal (antes de submeter scripts ao SLURM):**

```bash
module load Bio/FastQC/0.12.1
module load Bio/MultiQC/1.33
module load Bio/Cutadapt/5.2
```

**3. Configurar o acesso às bibliotecas R no RStudio Server:**

No início de cada sessão R (ou no topo do script), adicionar o caminho compartilhado:

```r
.libPaths(c("/opt/R/sharedLibs/4.3", .libPaths()))
```

**4. Instalar pacotes R (apenas na primeira vez):**

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("dada2", "phyloseq"))
install.packages(c("ggplot2", "vegan", "pheatmap", "reshape2", "gridExtra"))
```

### Parte 03 (Integração + MaAsLin3)

A Parte 03 **não** precisa de links simbólicos nem de módulos de terminal. Ela utiliza apenas R e depende dos **resultados já gerados** pelas Partes 01 e 02.

**Pré-requisitos:**

1. ✅ Ter concluído os pipelines DADA2 de ambas as partes (os arquivos `seqtab_nochim.rds` e `taxa.rds` devem existir em `parte_01/results/dada2/` e `parte_02/results/dada2/`).
2. ✅ Ter os metadados corretos em `parte_01/metadata.csv` e `parte_02/metadata.csv`.

**Pacote adicional (instalar apenas na primeira vez):**

```r
BiocManager::install("maaslin3")
```

**Como executar:**

```bash
cd parte_03
Rscript scripts/maaslin_analysis.R
```

> ⚠️ **Atenção:** Antes de executar, edite os caminhos genéricos dentro do script (`caminho/para/...`) para apontar para os arquivos corretos no seu sistema.

---

## 🛠️ Requisitos de Software (Resumo)

| Recurso | Versão | Usado em |
|---------|--------|----------|
| **R** | ≥ 4.3 | Todas as partes |
| `dada2` | Bioconductor | Partes 01, 02 e 03 |
| `phyloseq` | Bioconductor | Partes 01, 02 e 03 |
| `maaslin3` | Bioconductor | Parte 03 |
| `vegan`, `ggplot2`, `pheatmap` | CRAN | Partes 01 e 02 |
| `dplyr` | CRAN | Parte 03 |
| **FastQC** | 0.12.1 | Partes 01 e 02 |
| **MultiQC** | 1.33 | Partes 01 e 02 |
| **Cutadapt** | 5.2 | Partes 01 e 02 |

---

## 🚀 Como Começar

1. Clone o repositório no seu diretório do servidor:
   ```bash
   git clone https://github.com/cavalheiromf/LZT0693.git
   cd LZT0693
   ```
2. Siga a seção **Configuração Específica por Módulo** acima para preparar o ambiente da parte desejada.
3. Abra o `TUTORIAL.md` da parte correspondente e siga o passo a passo.

---

## 📜 Histórico de Modificações Recentes (Changelog)

Para fins de acompanhamento orgânico da evolução do projeto (especialmente para a Parte 03 e adaptações ao servidor):

* **Raiz do Projeto:**
  * Inclusão de instruções detalhadas no README sobre links simbólicos e carregamento de módulos (FastQC, MultiQC, Cutadapt).
  * O arquivo `.gitignore` foi corrigido para que todos os códigos fonte (`scripts/`) subam para o git, ignorando apenas saídas compiladas (`*.html`).
* **Módulo 3 (Integração e MaAsLin3):**
  * O tutorial foi totalmente reestruturado com sumário, diagramas e foco didático.
  * Adicionada a covariável `reads` (profundidade de sequenciamento) e a etapa de agrupamento taxonômico (`tax_glom`) para robustez estatística, seguindo boas práticas.
  * Criadas seções de Troubleshooting (Resolução de Problemas), Glossário e um exemplo concreto de interpretação dos resultados do MaAsLin3.
* **Módulos 1 e 2 (DADA2):**
  * A etapa de Aprendizado do Modelo de Erro (`learnErrors()`) foi simplificada. Como o servidor usa o R 4.3 (sem a função mais recente `makeBinnedQualErrfun`), o pipeline foi adaptado de forma didática para usar a função padrão do DADA2, conforme recomendação oficial para dados NovaSeq nessas condições.
