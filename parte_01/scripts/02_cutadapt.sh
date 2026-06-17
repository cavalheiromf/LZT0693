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

# Criar diretório de logs (se não existir)
mkdir -p logs

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
    "S813_16Sv3v4_01"
    "S813_16Sv3v4_02"
    "S813_16Sv3v4_03"
    "S813_16Sv3v4_04"
    "S813_16Sv3v4_05"
    "S813_16Sv3v4_06"
    "S813_16Sv3v4_07"
    "S813_16Sv3v4_NN"
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
