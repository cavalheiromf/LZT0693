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
