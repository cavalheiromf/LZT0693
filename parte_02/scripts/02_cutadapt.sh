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
