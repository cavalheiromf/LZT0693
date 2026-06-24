# Tutorial: Análise de Amplicon 16S rRNA (V3-V4) com DADA2

> **Projeto:** LZT0693 — Iniciação Científica em Biotecnologia  
> **Região alvo:** V3-V4 do gene 16S rRNA  
> **Plataforma:** Illumina (paired-end)  
> **Amostras:** 7 biológicas + 1 controle negativo

---

## Sumário

1. [Pré-requisitos e Configuração](#1-pré-requisitos-e-configuração)
2. [Controle de Qualidade (FastQC / MultiQC)](#2-controle-de-qualidade-fastqc--multiqc)
3. [Remoção de Primers (Cutadapt)](#3-remoção-de-primers-cutadapt)
4. [Pipeline DADA2 no R](#4-pipeline-dada2-no-r)
   - 4.1 [Preparação do ambiente](#41-preparação-do-ambiente)
   - 4.2 [Inspeção dos perfis de qualidade](#42-inspeção-dos-perfis-de-qualidade)
   - 4.3 [Filtragem e trimming](#43-filtragem-e-trimming)
   - 4.4 [Aprendizado do modelo de erro (Quality Binned)](#44-aprendizado-do-modelo-de-erro-quality-binned)
   - 4.5 [Inferência de ASVs (Denoising)](#45-inferência-de-asvs-denoising)
   - 4.6 [Merge dos pares e remoção de quimeras](#46-merge-dos-pares-e-remoção-de-quimeras)
   - 4.7 [Tracking das reads ao longo do pipeline](#47-tracking-das-reads-ao-longo-do-pipeline)
   - 4.8 [Atribuição taxonômica](#48-atribuição-taxonômica)
5. [Construção do objeto phyloseq](#5-construção-do-objeto-phyloseq)
6. [Diversidade Alfa](#6-diversidade-alfa)
7. [Diversidade Beta (Bray-Curtis)](#7-diversidade-beta-bray-curtis)
8. [Composição Taxonômica (Abundância Relativa)](#8-composição-taxonômica-abundância-relativa)

---

## 1. Pré-requisitos e Configuração

### Acesso ao RStudio Server

O RStudio Server está rodando no servidor e pode ser acessado diretamente pelo seu navegador web:

* **Endereço de Acesso:** [http://biotec02.esalq.usp.br:59087](http://biotec02.esalq.usp.br:59087)

> **💡 Dica para Acesso Remoto (Túnel SSH):**  
> Se você estiver fora da rede da USP (ou se o acesso direto estiver bloqueado pelo firewall), você pode criar um túnel SSH a partir do terminal do seu computador pessoal:
> ```bash
> ssh -L 59087:localhost:59087 seu_usuario@biotec02.esalq.usp.br
> ```
> E então acesse no seu navegador: **[http://localhost:59087](http://localhost:59087)**

---

### Software no servidor (via `module load`)

As ferramentas de linha de comando estão disponíveis no servidor via sistema de módulos.
Para carregar as versões mais recentes:

```bash
# Carregar módulos necessários
module load Bio/FastQC/0.12.1
module load Bio/MultiQC/1.33
module load Bio/Cutadapt/5.2
```

> **Dica:** Para verificar todas as versões disponíveis, use `module avail Bio/FastQC`, etc.
> Para listar os módulos carregados: `module list`

### Pacotes R

Os pacotes R devem ser instalados no RStudio ou em uma sessão R interativa:

```r
# Instalar BiocManager se necessário
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Pacotes do Bioconductor
BiocManager::install(c("dada2", "phyloseq"))

# Pacotes do CRAN
install.packages(c("ggplot2", "vegan", "pheatmap", "reshape2", "gridExtra"))
```

### Banco de dados taxonômico

O banco SILVA formatado para o DADA2 está disponível diretamente no servidor. Você criará um link simbólico para a pasta contendo os bancos (veja na seção abaixo), que contém os seguintes arquivos da versão v138.1:
- `silva_nr99_v138.1_train_set.fa.gz`
- `silva_species_assignment_v138.1.fa.gz`

### O Sistema de Gerenciamento de Filas (SLURM)

Como o processamento inicial de bioinformática (FastQC, MultiQC e Cutadapt) exige alto poder computacional, ele é executado no servidor utilizando o gerenciador de recursos **SLURM**. 

O SLURM gerencia o uso compartilhado do hardware (CPU, memória e armazenamento) entre todos os usuários, garantindo que o servidor não sofra sobrecargas.

#### Comandos essenciais do SLURM:

1. **Submeter um script de processamento**:
   ```bash
   sbatch scripts/meu_script.sh
   ```
2. **Acompanhar o status das tarefas em execução ou na fila**:
   ```bash
   # Listar todas as tarefas ativas no servidor
   squeue
   
   # Listar apenas as suas tarefas ativas
   squeue -u $USER
   ```
3. **Cancelar uma tarefa enviada por engano**:
   ```bash
   scancel <JOB_ID>
   ```
4. **Verificar o estado das partições (filas) do cluster**:
   ```bash
   sinfo
   ```

#### Estrutura do cabeçalho dos scripts SLURM:
Nossos scripts Bash (`scripts/*.sh`) contêm diretivas especiais `#SBATCH` no início. Elas informam ao gerenciador os recursos necessários:
* `#SBATCH --partition=short`: Envia a tarefa para a fila rápida `short` (limite de 1 dia de execução).
* `#SBATCH --cpus-per-task=4`: Solicita 4 núcleos de processamento paralelos para acelerar a execução (ex: threads no FastQC/Cutadapt).
* `#SBATCH --mem=8G`: Reserva 8 Gigabytes de memória RAM para que o programa não sofra interrupção por falta de memória.
* `#SBATCH --output=logs/...`: Define onde salvar as saídas de texto e relatórios gerados pelos programas.

> **💡 Informações e FAQ do Servidor:**  
> Para consultar limites de tempo das partições, configurações físicas do cluster e outras informações práticas do servidor utilizado em aula, acesse a página de FAQ oficial:  
> 🔗 **[FAQ do Servidor Biotec02](http://biotec02.esalq.usp.br:59080/faq.html)**

---

## Estrutura de diretórios esperada (antes da análise)

Arquivos obtidos ao clonar o repositório:

```
LZT0693/
├── slides/             # Slides das aulas de bioinformática
│   └── Aula_Do_FASTQ_a_Taxonomia_16S_Completa.pdf
├── metadata.csv        # Metadados das amostras
├── TUTORIAL.md         # Este tutorial
├── .gitignore
└── README.md
```

### Como montar a estrutura de trabalho?

Apenas dois links simbólicos precisam ser criados manualmente antes de iniciar o pipeline:

```bash
# 1. Dados brutos de sequenciamento (FASTQ)
ln -s /home/cursos/LCoutinho202604/data_shared/aula_03/data data

# 2. Banco de dados taxonômico SILVA
ln -s /home/cursos/LCoutinho202604/data_shared/aula_03/databases databases
```

O restante da estrutura é gerado automaticamente durante o pipeline:

- **`scripts/`** — crie os arquivos conforme as seções 2, 3 e 4 deste tutorial (ou cole os scripts fornecidos em aula)
- **`results/`** e subpastas — criados pelos comandos `mkdir -p` nos scripts Bash e `dir.create()` no R
- **`logs/`** — criado pelos scripts SLURM ao serem submetidos

---

## 2. Controle de Qualidade (FastQC / MultiQC)

Crie o arquivo `scripts/01_fastqc.sh` e submeta via SLURM:

```bash
#!/bin/bash
#SBATCH --job-name=fastqc_16S
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fastqc_%j.out
#SBATCH --error=logs/fastqc_%j.err

# ==============================================================
# Controle de qualidade com FastQC e MultiQC
# ==============================================================

# Carregar módulos
module load Bio/FastQC/0.12.1
module load Bio/MultiQC/1.33

# Diretórios
DATADIR="data"
OUTDIR="results/fastqc"
MULTIQCDIR="results/multiqc"

# Criar diretórios de saída
mkdir -p ${OUTDIR}

# Lista de amostras
SAMPLES=(
    "S813_A1"
    "S813_A2"
    "S813_A3"
    "S813_A4"
    "S813_A5"
    "S813_A6"
    "S813_A7"
    "S813_A8"
    "S813_A9"
    "S813_A10"
    "S813_A11"
    "S813_A12"
    "S813_B1"
    "S813_B2"
    "S813_B3"
    "S813_B4"
    "S813_NN_V3V4"
)

# Executar FastQC para cada amostra
for SAMPLE in "${SAMPLES[@]}"; do
    echo ">> FastQC: ${SAMPLE}"
    mkdir -p ${OUTDIR}/${SAMPLE}
    fastqc \
        ${DATADIR}/${SAMPLE}/*.fastq.gz \
        --outdir ${OUTDIR}/${SAMPLE} \
        --threads ${SLURM_CPUS_PER_TASK:-4}
done

# Gerar relatório consolidado com MultiQC
echo ">> MultiQC: gerando relatório consolidado"
mkdir -p ${MULTIQCDIR}
multiqc ${OUTDIR} \
    --outdir ${MULTIQCDIR} \
    --force

echo ">> FastQC/MultiQC finalizado."
```

Submeter ao SLURM:

```bash
mkdir -p logs
sbatch scripts/01_fastqc.sh
```

**O que verificar no relatório:**

- Qualidade média por posição (Phred ≥ 20 é aceitável; ≥ 30 é bom)
- Conteúdo de adaptadores
- Distribuição de tamanho de reads
- Nível de duplicação
- O controle negativo (NN) deve apresentar muito poucos reads

---

## 3. Remoção de Primers (Cutadapt)

As sequências de primers V3-V4 mais comuns são:

| Primer | Sequência (5'→3') | Tamanho |
|---|---|---|
| 341F (forward) | `CCTACGGGNGGCWGCAG` | 17 nt |
| 805R (reverse) | `GACTACHVGGGTATCTAATCC` | 21 nt |

> **⚠️ Importante:** Confirme as sequências de primers utilizadas no seu experimento antes de executar esta etapa. As sequências acima são para os primers 341F/805R comumente usados para V3-V4, mas podem variar.

Crie o arquivo `scripts/02_cutadapt.sh` e submeta via SLURM:

```bash
#!/bin/bash
#SBATCH --job-name=cutadapt_16S
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cutadapt_%j.out
#SBATCH --error=logs/cutadapt_%j.err

# ==============================================================
# Remoção de primers com Cutadapt
# ==============================================================

# Carregar módulo
module load Bio/Cutadapt/5.2

# Sequências dos primers (CONFIRMAR antes de executar!)
FWD="CCTACGGGNGGCWGCAG"       # 341F
REV="GACTACHVGGGTATCTAATCC"   # 805R

# Complemento reverso dos primers (para paired-end)
FWD_RC="CTGCWGCCNCCGTAGG"
REV_RC="GGATTAGATACCCBDGTAGTC"

# Diretórios
DATADIR="data"
OUTDIR="results/cutadapt"
LOGDIR="results/cutadapt/logs"

mkdir -p ${OUTDIR} ${LOGDIR}

# Lista de amostras
SAMPLES=(
    "S813_A1"
    "S813_A2"
    "S813_A3"
    "S813_A4"
    "S813_A5"
    "S813_A6"
    "S813_A7"
    "S813_A8"
    "S813_A9"
    "S813_A10"
    "S813_A11"
    "S813_A12"
    "S813_B1"
    "S813_B2"
    "S813_B3"
    "S813_B4"
    "S813_NN_V3V4"
)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "=============================================="
    echo ">> Cutadapt: ${SAMPLE}"
    echo "=============================================="

    # Identificar os arquivos R1 e R2
    R1=$(ls ${DATADIR}/${SAMPLE}/*_R1_001.fastq.gz)
    R2=$(ls ${DATADIR}/${SAMPLE}/*_R2_001.fastq.gz)

    # Saída
    OUT_R1="${OUTDIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    OUT_R2="${OUTDIR}/${SAMPLE}_R2_trimmed.fastq.gz"

    cutadapt \
        -g "${FWD}" \
        -a "${REV_RC}" \
        -G "${REV}" \
        -A "${FWD_RC}" \
        --discard-untrimmed \
        --minimum-length 100 \
        -j ${SLURM_CPUS_PER_TASK:-4} \
        -o "${OUT_R1}" \
        -p "${OUT_R2}" \
        "${R1}" "${R2}" \
        > "${LOGDIR}/${SAMPLE}_cutadapt.log" 2>&1

    echo "   Reads processados. Log: ${LOGDIR}/${SAMPLE}_cutadapt.log"
done

echo ""
echo ">> Cutadapt finalizado para todas as amostras."
echo ">> Verifique os logs em: ${LOGDIR}/"
```

Submeter ao SLURM:

```bash
sbatch scripts/02_cutadapt.sh
```

**O que verificar nos logs do Cutadapt:**

- % de reads com primers encontrados (deve ser alta, >80%)
- Reads descartados (sem primer) — devem ser poucos
- O controle negativo pode ter pouquíssimos reads com primers

---

## 4. Pipeline DADA2 no R

A partir daqui, todo o processamento é feito em **R**. Execute interativamente no RStudio.

### 4.1 Preparação do ambiente

```r
# ==============================================================
# Pipeline DADA2 — Análise de amplicon 16S rRNA (V3-V4)
# Projeto: LZT0693 — Iniciação Científica em Biotecnologia
# ==============================================================

# Para garantir o acesso às bibliotecas de bioinformática pré-instaladas no servidor,
# adicione o caminho compartilhado ao início dos seus locais de busca (.libPaths)
.libPaths(c("/opt/R/sharedLibs/4.3", .libPaths()))

library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(pheatmap)

# Diretório base do projeto (ajustar conforme necessário)
base_dir <- "."   # ou caminho absoluto para S813_aulaBiotec26/

# Diretório com reads trimados pelo Cutadapt
cut_dir <- file.path(base_dir, "results", "cutadapt")

# Diretório de saída DADA2
dada2_dir <- file.path(base_dir, "results", "dada2")
dir.create(dada2_dir, recursive = TRUE, showWarnings = FALSE)

# Listar arquivos R1 e R2 (trimados)
fnFs <- sort(list.files(cut_dir, pattern = "_R1_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(cut_dir, pattern = "_R2_trimmed.fastq.gz", full.names = TRUE))

# Extrair nomes das amostras
sample_names <- gsub("_R1_trimmed.fastq.gz", "", basename(fnFs))
cat("Amostras encontradas:", length(sample_names), "\n")
print(sample_names)
```

### 4.2 Inspeção dos perfis de qualidade

```r
# Perfis de qualidade — Forward (R1)
p_fwd <- plotQualityProfile(fnFs, aggregate = TRUE)
ggsave(file.path(dada2_dir, "quality_profile_R1.png"), p_fwd,
       width = 10, height = 6, dpi = 150)

# Perfis de qualidade — Reverse (R2)
p_rev <- plotQualityProfile(fnRs, aggregate = TRUE)
ggsave(file.path(dada2_dir, "quality_profile_R2.png"), p_rev,
       width = 10, height = 6, dpi = 150)
```

> **Dica:** Use estes gráficos para decidir os pontos de truncamento (`truncLen`).
> Tipicamente, para dados V3-V4 Illumina 2×300:
> - R1: truncar em ~280 nt (qualidade geralmente boa)
> - R2: truncar em ~200-250 nt (qualidade cai mais cedo)
>
> Os reads precisam ter overlap suficiente para o merge (~20 nt mínimo).
> Amplicon V3-V4 ≈ 460 bp → truncF + truncR deve cobrir >460 + 20 = 480 bp.

> **⚠️ Nota sobre `truncLen` — NovaSeq vs MiSeq:**
>
> A estratégia de truncamento depende da **plataforma de sequenciamento**:
>
> | Plataforma | Quality scores | Estratégia de `truncLen` |
> |---|---|---|
> | **MiSeq** | Phred scores contínuos (tradicionais) | Definir `truncLen` com base nos perfis de qualidade (ex: `c(280, 200)`). A qualidade cai progressivamente e o truncamento é essencial. |
> | **NovaSeq** (e NextSeq/NovaSeq X) | Quality scores **binned** (apenas uns poucos valores discretos: ~2, 12, 23, 37) | Pode-se usar `truncLen = c(0, 0)` (sem truncamento por posição) e confiar no filtro `truncQ = 8` para remover caudas de baixa qualidade. |
>
> **Por que a diferença?** O NovaSeq usa **quality score binning** (RTA3),
> o que significa que os quality scores são arredondados para poucos valores
> discretos. Isso torna os perfis de qualidade "em escada" e dificulta a
> escolha visual de um ponto de truncamento. Para esses dados, a abordagem
> recomendada é:
>
> ```r
> # Para dados NovaSeq / binned quality:
> truncLen = c(0, 0)   # Não truncar por posição
> truncQ = 8           # Truncar quando Phred cair abaixo de 8
> ```
>
> **Verifique seus perfis de qualidade:** se os scores aparecem em "degraus"
> (ex: blocos em Q37, Q23, Q12), seus dados provavelmente são do NovaSeq
> e você deve usar a estratégia com `truncQ` em vez de `truncLen` fixo.

### 4.3 Filtragem e trimming

```r
# ==============================================================
# Filtragem e trimming
# ==============================================================

# Diretório de saída para reads filtradas
filt_dir <- file.path(dada2_dir, "filtered")
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)

# Caminhos dos arquivos filtrados
filtFs <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# ---------------------------------------------------------------
# Parâmetros de truncamento
# ---------------------------------------------------------------
# OPÇÃO A — MiSeq (Phred scores contínuos):
#   Ajustar com base nos perfis de qualidade (seção 4.2)
#truncF <- 280
#truncR <- 200
#truncQ_val <- 2

# OPÇÃO B — NovaSeq (quality scores binned):
#   Descomentar as linhas abaixo e comentar a Opção A
truncF <- 0
truncR <- 0
truncQ_val <- 8
# ---------------------------------------------------------------

# Executar a filtragem
filt_out <- filterAndTrim(
    fnFs, filtFs,
    fnRs, filtRs,
    truncLen = c(truncF, truncR),
    maxN = 0,
    maxEE = c(2, 2),
    truncQ = truncQ_val,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
)

# Estatísticas de filtragem
filt_stats <- as.data.frame(filt_out)
filt_stats$pct_passed <- round(filt_stats$reads.out / filt_stats$reads.in * 100, 1)

cat("\n>> Estatísticas de filtragem:\n")
print(filt_stats)

# Salvar
write.csv(filt_stats, file.path(filt_dir, "filter_stats.csv"), row.names = TRUE)
```

### 4.4 Aprendizado do modelo de erro

O DADA2 precisa estimar as taxas de erro de sequenciamento. Utilizaremos a função padrão `learnErrors()`.

> **⚠️ Nota sobre R 4.3 e dados NovaSeq (Quality Binned):**
> O NovaSeq gera *quality scores binned* (agrupados em degraus). Versões mais recentes do DADA2 possuem uma função específica para isso (`makeBinnedQualErrfun`). Porém, como estamos utilizando um ambiente com **R 4.3** no servidor e uma versão mais antiga do pacote, essa função não está disponível.
> 
> A solução recomendada pelo próprio desenvolvedor do DADA2 para essas versões antigas é **utilizar a função padrão `learnErrors`**. O modelo padrão (loess) tentará se ajustar aos degraus. Os pontos nos gráficos de erro parecerão estranhos (agrupados verticalmente em valores discretos de Phred), mas a linha preta de regressão ainda assim estimará as taxas de erro corretamente.

```r
# ==============================================================
# Aprendizado do modelo de erro
# ==============================================================

# Aprender o modelo de erro — Forward
cat(">> Aprendendo modelo de erro (Forward)...\n")
errF <- learnErrors(filtFs, multithread = TRUE)

# Aprender o modelo de erro — Reverse
cat(">> Aprendendo modelo de erro (Reverse)...\n")
errR <- learnErrors(filtRs, multithread = TRUE)

# Visualizar os modelos de erro
p_errF <- plotErrors(errF, nominalQ = TRUE)
ggsave(file.path(dada2_dir, "error_model_Forward.png"), p_errF,
       width = 10, height = 8, dpi = 150)

p_errR <- plotErrors(errR, nominalQ = TRUE)
ggsave(file.path(dada2_dir, "error_model_Reverse.png"), p_errR,
       width = 10, height = 8, dpi = 150)
```

> **Verificação:** Os pontos pretos nos gráficos de erro devem acompanhar
> razoavelmente as linhas pretas (taxas de erro estimadas). Se houver grandes
> discrepâncias, considere ajustar os parâmetros de filtragem.

### 4.5 Inferência de ASVs (Denoising)

```r
# ==============================================================
# Inferência de ASVs (Denoising)
# ==============================================================

cat(">> Inferindo ASVs (Forward)...\n")
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)

cat(">> Inferindo ASVs (Reverse)...\n")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Resumo
cat("\nASVs por amostra (Forward):\n")
print(sapply(dadaFs, function(x) length(x$denoised)))
```

### 4.6 Merge dos pares e remoção de quimeras

```r
# ==============================================================
# Merge das reads paired-end
# ==============================================================

cat(">> Realizando merge dos pares R1/R2...\n")
merged <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Construir tabela de sequências
seqtab <- makeSequenceTable(merged)
cat("\nDimensões da tabela de sequências:", dim(seqtab), "\n")

# Distribuição de tamanhos dos amplicons
cat("Distribuição de tamanhos:\n")
print(table(nchar(getSequences(seqtab))))

# ==============================================================
# Remoção de quimeras
# ==============================================================

cat(">> Removendo quimeras...\n")
seqtab_nochim <- removeBimeraDenovo(seqtab,
                                      method = "consensus",
                                      multithread = TRUE,
                                      verbose = TRUE)

cat("\nASVs após remoção de quimeras:", ncol(seqtab_nochim), "\n")
cat("Fração de reads retidos:",
    round(sum(seqtab_nochim) / sum(seqtab) * 100, 1), "%\n")
```

### 4.7 Tracking das reads ao longo do pipeline

```r
# ==============================================================
# Rastreamento de reads por etapa
# ==============================================================

getN <- function(x) sum(getUniques(x))

track <- data.frame(
    input      = filt_out[, 1],
    filtered   = filt_out[, 2],
    denoisedF  = sapply(dadaFs, getN),
    denoisedR  = sapply(dadaRs, getN),
    merged     = sapply(merged, getN),
    nonchim    = rowSums(seqtab_nochim)
)
rownames(track) <- sample_names

cat("\n>> Tracking de reads:\n")
print(track)

# Salvar
write.csv(track, file.path(dada2_dir, "read_tracking.csv"))
```

### 4.8 Atribuição taxonômica

```r
# ==============================================================
# Atribuição taxonômica com SILVA
# ==============================================================

# Caminhos para o banco SILVA (EDITAR conforme necessário)
silva_train <- file.path(base_dir, "databases", "silva_nr99_v138.1_train_set.fa.gz")
silva_species <- file.path(base_dir, "databases", "silva_species_assignment_v138.1.fa.gz")

cat(">> Atribuindo taxonomia (pode demorar alguns minutos)...\n")
taxa <- assignTaxonomy(seqtab_nochim, silva_train, multithread = TRUE)

cat(">> Adicionando atribuição a nível de espécie...\n")
taxa <- addSpecies(taxa, silva_species)

# Inspecionar
cat("\nPrimeiras linhas da tabela taxonômica:\n")
taxa_print <- taxa
rownames(taxa_print) <- NULL  # Para visualização mais limpa
print(head(taxa_print))

# Salvar
saveRDS(seqtab_nochim, file.path(dada2_dir, "seqtab_nochim.rds"))
saveRDS(taxa, file.path(dada2_dir, "taxa.rds"))
write.csv(taxa, file.path(dada2_dir, "taxonomy_table.csv"))
```

---

## 5. Construção do objeto phyloseq

```r
# ==============================================================
# Construir objeto phyloseq
# ==============================================================

library(phyloseq)

# Carregar metadados
metadata <- read.csv(file.path(base_dir, "metadata.csv"), row.names = 1)

# Verificar que os nomes das amostras batem
cat("Amostras nos dados:", rownames(seqtab_nochim), "\n")
cat("Amostras nos metadados:", rownames(metadata), "\n")

# Criar o objeto phyloseq
ps <- phyloseq(
    otu_table(seqtab_nochim, taxa_are_rows = FALSE),
    sample_data(metadata),
    tax_table(taxa)
)

# Renomear ASVs para nomes mais legíveis
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

cat("\n>> Objeto phyloseq criado:\n")
print(ps)

# (Opcional) Remover o controle negativo antes das análises de diversidade
# ps <- prune_samples(sample_names(ps) != "S813_NN_V3V4", ps)

# Salvar
saveRDS(ps, file.path(dada2_dir, "phyloseq_object.rds"))
```

---

## 6. Diversidade Alfa

```r
# ==============================================================
# Diversidade Alfa
# ==============================================================

# Remover controle negativo para análises ecológicas
ps_bio <- prune_samples(sample_names(ps) != "S813_NN_V3V4", ps)

# Remover ASVs com zero reads após remoção do controle
ps_bio <- prune_taxa(taxa_sums(ps_bio) > 0, ps_bio)

# Calcular índices de diversidade alfa
alpha_div <- estimate_richness(ps_bio, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
alpha_div$sample_id <- rownames(alpha_div)

cat("\n>> Diversidade Alfa:\n")
print(alpha_div)

# Salvar tabela
write.csv(alpha_div, file.path(dada2_dir, "alpha_diversity.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# Gráfico de diversidade alfa
# ------------------------------------------------------------------

# Adicionar metadados
alpha_div_meta <- merge(alpha_div, data.frame(sample_data(ps_bio)), by = "row.names")

# Gráfico: Diversidade Alfa por Grupo
p_obs <- plot_richness(ps_bio, x = "group", measures = c("Observed", "Shannon", "Simpson"), color = "group") +
    geom_boxplot(alpha = 0.3) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(title = "Diversidade Alfa por Grupo — 16S V3-V4")

ggsave(file.path(dada2_dir, "alpha_diversity_plot.png"), p_obs,
       width = 12, height = 6, dpi = 150)

cat(">> Gráficos de diversidade alfa salvos.\n")
```

> **Interpretação dos índices:**
>
> | Índice | O que mede | Valores |
> |--------|------------|---------|
> | **Observed** | Número de ASVs (riqueza) | Inteiro positivo |
> | **Chao1** | Riqueza estimada (inclui raras) | ≥ Observed |
> | **Shannon** | Diversidade (riqueza + equitabilidade) | 0 → ∞ (maior = mais diverso) |
> | **Simpson** | Dominância (1 - probabilidade de sortear 2 iguais) | 0 → 1 (maior = mais diverso) |

---

## 7. Diversidade Beta (Bray-Curtis)

```r
# ==============================================================
# Diversidade Beta — Bray-Curtis
# ==============================================================

library(vegan)
library(pheatmap)

# Transformar para proporções (abundância relativa)
ps_rel <- transform_sample_counts(ps_bio, function(x) x / sum(x))

# ------------------------------------------------------------------
# 7.1 Calcular distância de Bray-Curtis
# ------------------------------------------------------------------

bray_dist <- phyloseq::distance(ps_rel, method = "bray")

cat("\n>> Matriz de distância de Bray-Curtis:\n")
print(as.matrix(bray_dist))

# Salvar a matriz
write.csv(as.matrix(bray_dist),
          file.path(dada2_dir, "bray_curtis_matrix.csv"))

# ------------------------------------------------------------------
# 7.2 PCoA (Principal Coordinates Analysis)
# ------------------------------------------------------------------

pcoa_res <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)

# Gráfico PCoA colorido por grupo
p_pcoa <- plot_ordination(ps_rel, pcoa_res, type = "samples", color = "group") +
    geom_point(size = 4, alpha = 0.8) +
    theme_bw() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12)
    ) +
    labs(title = "PCoA — Bray-Curtis — 16S V3-V4") +
    stat_ellipse(level = 0.95, linetype = 2)

ggsave(file.path(dada2_dir, "pcoa_bray_curtis.png"), p_pcoa,
       width = 8, height = 6, dpi = 150)

# ------------------------------------------------------------------
# 7.3 Heatmap da matriz de Bray-Curtis
# ------------------------------------------------------------------

bray_mat <- as.matrix(bray_dist)

pheatmap(bray_mat,
         clustering_distance_rows = bray_dist,
         clustering_distance_cols = bray_dist,
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
         main = "Bray-Curtis — Heatmap",
         fontsize = 10,
         filename = file.path(dada2_dir, "bray_curtis_heatmap.png"),
         width = 8, height = 7)

# ------------------------------------------------------------------
# 7.4 PERMANOVA (teste estatístico)
# ------------------------------------------------------------------

# Extrair metadados como data.frame
meta_df <- data.frame(sample_data(ps_bio))

# PERMANOVA — testar se os grupos diferem
set.seed(42)
permanova_res <- adonis2(bray_dist ~ group,
                          data = meta_df,
                          permutations = 999)
cat("\n>> Resultado PERMANOVA:\n")
print(permanova_res)
write.csv(as.data.frame(permanova_res),
          file.path(dada2_dir, "permanova_bray_curtis.csv"))

cat("\n>> Análise de diversidade beta concluída.\n")
cat(">> Resultados salvos em:", dada2_dir, "\n")
```

> **Interpretação:**
>
> - **Bray-Curtis:** Mede dissimilaridade baseada em abundância.
>   Valores de 0 (idênticas) a 1 (completamente diferentes).
> - **PCoA:** Reduz a matriz de distâncias a eixos principais.
>   Amostras próximas no gráfico são mais similares.
> - **PERMANOVA:** Testa se a composição difere significativamente entre grupos
>   (p < 0.05 indica diferença significativa).

---

## 8. Composição Taxonômica (Abundância Relativa)

Além de índices ecológicos gerais, é fundamental visualizar a composição real de microrganismos em cada amostra. Com o `phyloseq`, podemos facilmente colapsar a taxonomia por Filo e plotar um gráfico de barras acumuladas.

```r
# ==============================================================
# Composição Taxonômica por Filo
# ==============================================================

# Agrupar/aglomerar por nível de Filo
ps_phylum <- tax_glom(ps_bio, "Phylum")

# Converter abundância absoluta para abundância relativa
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

# Plotar gráfico de barras acumuladas de abundância relativa
p_bar <- plot_bar(ps_phylum_rel, fill = "Phylum") +
    geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(size = 14, face = "bold")
    ) +
    labs(title = "Composição Taxonômica por Filo (Abundância Relativa)",
         x = "Amostra",
         y = "Abundância Relativa")

ggsave(file.path(dada2_dir, "taxonomic_composition_phylum.png"), p_bar,
       width = 10, height = 6, dpi = 150)
```

> **Interpretação:**
> O gráfico mostra a abundância relativa acumulada de filos bacterianos em cada amostra (ou grupo), facilitando a comparação direta do perfil de táxons predominantes.

---

## Estrutura de diretórios esperada (após a análise)

Esta é a estrutura completa que você terá localmente após executar todas as etapas. Os itens marcados com `*` **não estão no repositório** — são criados ou linkados localmente.

```
LZT0693/
│
│  # --- Rastreados pelo git ---
├── slides/
│   └── Aula_Do_FASTQ_a_Taxonomia_16S_Completa.pdf
├── metadata.csv
├── TUTORIAL.md
├── .gitignore
├── README.md
│
│  # --- Criados/linkados localmente (não rastreados pelo git) ---
├── data/ *              → symlink para dados brutos no servidor
│   ├── S813_A1/
│   │   ├── S813_A1_S1_L001_R1_001.fastq.gz
│   │   └── S813_A1_S1_L001_R2_001.fastq.gz
│   ├── S813_A2/ ... S813_B4/
│   └── S813_NN_V3V4/
├── databases/ *        → symlink para bancos SILVA no servidor
│   ├── silva_nr99_v138.1_train_set.fa.gz
│   └── silva_species_assignment_v138.1.fa.gz
├── scripts/ *          # Scripts criados conforme o tutorial
│   ├── 01_fastqc.sh
│   ├── 02_cutadapt.sh
│   └── microbiome_tutorial.R
├── logs/ *             # Gerado pelo SLURM
│   ├── fastqc_*.out / .err
│   └── cutadapt_*.out / .err
└── results/ *          # Gerado pelo pipeline
    ├── fastqc/                              # Etapa 2
    │   ├── S813_A1/
    │   │   ├── *_R1_001_fastqc.html
    │   │   ├── *_R1_001_fastqc.zip
    │   │   ├── *_R2_001_fastqc.html
    │   │   └── *_R2_001_fastqc.zip
    │   └── ... (por amostra)
    ├── multiqc/                             # Etapa 2
    │   ├── multiqc_report.html
    │   └── multiqc_data/
    ├── cutadapt/                            # Etapa 3
    │   ├── S813_16Sv3v4_01_R1_trimmed.fastq.gz
    │   ├── S813_16Sv3v4_01_R2_trimmed.fastq.gz
    │   ├── ... (por amostra)
    │   └── logs/
    │       └── S813_16Sv3v4_01_cutadapt.log ...
    └── dada2/                               # Etapas 4–8
        ├── filtered/
        │   ├── S813_16Sv3v4_01_F_filt.fastq.gz
        │   ├── S813_16Sv3v4_01_R_filt.fastq.gz
        │   ├── ... (por amostra)
        │   └── filter_stats.csv
        ├── quality_profile_R1.png
        ├── quality_profile_R2.png
        ├── error_model_Forward.png
        ├── error_model_Reverse.png
        ├── read_tracking.csv
        ├── seqtab_nochim.rds
        ├── taxa.rds
        ├── taxonomy_table.csv
        ├── phyloseq_object.rds
        ├── alpha_diversity.csv
        ├── alpha_diversity_plot.png
        ├── bray_curtis_matrix.csv
        ├── bray_curtis_heatmap.png
        ├── pcoa_bray_curtis.png
        └── permanova_bray_curtis.csv
```

---

## Próximos passos (opcionais)

- [ ] Ajustar `truncLen` com base nos perfis de qualidade e na plataforma (MiSeq vs NovaSeq)
- [ ] Confirmar os primers utilizados no experimento
- [ ] Executar PERMANOVA para testar diferenças entre grupos
- [ ] Explorar composição taxonômica (barplots por Filo/Gênero)

---

## Referências

1. **DADA2:** Callahan, B.J. et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583. [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
2. **SILVA:** Quast, C. et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*, 41, D590–D596. [doi:10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
3. **phyloseq:** McMurdie, P.J. & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis. *PLoS ONE*, 8(4), e61217. [doi:10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)
4. **Cutadapt:** Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, 17(1), 10–12. [doi:10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)
5. **vegan:** Oksanen, J. et al. (2022). vegan: Community Ecology Package. [CRAN](https://CRAN.R-project.org/package=vegan)
