# LZT0693 — Análise de Amplicon 16S rRNA (V3-V4) | Iniciação Científica em Biotecnologia

## Descrição

Este repositório contém os dados brutos e o pipeline de análise de sequenciamento de amplicon **16S rRNA (regiões V3-V4)** gerados no contexto da disciplina LZT0693 - Iniciação Científica em Biotecnologia.

O sequenciamento foi realizado em plataforma **Illumina** (paired-end), resultando em 8 bibliotecas: 7 amostras biológicas e 1 controle negativo.

## Como Clonar / Baixar este Repositório

Para clonar este repositório no terminal do servidor, execute:

```bash
git clone https://github.com/cavalheiromf/LZT0693.git
cd LZT0693
```

## Estrutura do repositório

```
LZT0693/
├── data/                          # Dados brutos (FASTQ)
│   ├── S813_16Sv3v4_01/           # Amostra 01 (R1 + R2)
│   ├── ...                        # Demais amostras
│   └── S813_16Sv3v4_NN/           # Controle negativo
├── scripts/                       # Scripts da análise (.sh e .R)
│   ├── 01_fastqc.sh               # Controle de qualidade
│   ├── 02_cutadapt.sh             # Remoção de primers
│   └── microbiome_tutorial.R      # Pipeline DADA2 e ecologia em R
├── results/                       # Resultados das análises
│   ├── fastqc/                    # Relatórios FastQC por amostra
│   └── multiqc/                   # Relatório consolidado MultiQC
├── metadata.csv                   # Metadados das amostras (editar!)
├── TUTORIAL.md                    # Tutorial completo do pipeline
└── README.md                      # Este arquivo
```

## Amostras

| ID da Amostra | Index Illumina | Tipo |
|---|---|---|
| S813_16Sv3v4_01 | S83 | Amostra biológica |
| S813_16Sv3v4_02 | S84 | Amostra biológica |
| S813_16Sv3v4_03 | S85 | Amostra biológica |
| S813_16Sv3v4_04 | S86 | Amostra biológica |
| S813_16Sv3v4_05 | S87 | Amostra biológica |
| S813_16Sv3v4_06 | S88 | Amostra biológica |
| S813_16Sv3v4_07 | S89 | Amostra biológica |
| S813_16Sv3v4_NN | S90 | Controle negativo |

## Pipeline de Análise

O pipeline completo está descrito em detalhes no arquivo [`TUTORIAL.md`](TUTORIAL.md) e inclui:

1. **Controle de qualidade** — FastQC + MultiQC
2. **Remoção de primers** — Cutadapt
3. **Denoising e inferência de ASVs** — DADA2 (com modelo de erro *quality binned*: bins 12, 24, 40)
4. **Atribuição taxonômica** — DADA2 + banco SILVA
5. **Diversidade alfa** — Índices de riqueza e diversidade (Shannon, Simpson, Observed, Chao1)
6. **Diversidade beta** — Dissimilaridade de Bray-Curtis + PCoA

## Pré-requisitos

- **R** ≥ 4.3 com os pacotes: `dada2`, `phyloseq`, `vegan`, `ggplot2`, `pheatmap`
- **Cutadapt** 5.2 (via `module load Bio/Cutadapt/5.2`)
- **FastQC** 0.12.1 (via `module load Bio/FastQC/0.12.1`)
- **MultiQC** 1.33 (via `module load Bio/MultiQC/1.33`)
- Acesso ao cluster SLURM para submissão dos scripts bash

## Como começar

1. Edite o arquivo `metadata.csv` com as informações reais das amostras
2. Siga o tutorial passo a passo em `TUTORIAL.md`

## Referências

- Callahan, B.J. et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583.
- Quast, C. et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*, 41, D590–D596.
- McMurdie, P.J. & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE*, 8(4), e61217.

## Licença

Material didático para uso interno no curso de Iniciação Científica em Biotecnologia - LZT0693.
