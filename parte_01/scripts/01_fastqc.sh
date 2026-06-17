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
    "S813_16Sv3v4_01"
    "S813_16Sv3v4_02"
    "S813_16Sv3v4_03"
    "S813_16Sv3v4_04"
    "S813_16Sv3v4_05"
    "S813_16Sv3v4_06"
    "S813_16Sv3v4_07"
    "S813_16Sv3v4_NN"
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
