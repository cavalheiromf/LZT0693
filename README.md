# LZT0693 — Análise de Amplicon 16S rRNA (V3-V4) | Iniciação Científica em Biotecnologia

## Descrição

Este repositório contém o pipeline de análise de sequenciamento de amplicon **16S rRNA (regiões V3-V4)** desenvolvido na disciplina LZT0693 - Iniciação Científica em Biotecnologia.

O sequenciamento foi realizado em plataforma **Illumina** (paired-end), resultando em 8 bibliotecas: 7 amostras biológicas (líquen, musgo e mel) e 1 controle negativo.

> **Nota:** Os dados brutos (FASTQ) **não estão incluídos** neste repositório. Eles são acessados por meio de um link simbólico (`data/`) apontando para o diretório compartilhado do servidor.

## Como Clonar / Baixar este Repositório

Para clonar este repositório no terminal do servidor, execute:

```bash
git clone https://gitlab.com/cavalheiromf/LZT0693.git
cd LZT0693
```

Após clonar, recriar os links simbólicos necessários:

```bash
ln -s /home/cursos/LCoutinho202604/data_shared/aula_02/data     data
ln -s /home/cursos/LCoutinho202604/data_shared/aula_02/databases databases
```

## Estrutura do Repositório

Arquivos rastreados pelo git (o que você obtém ao clonar):

```
LZT0693/
├── slides/             # Slides das aulas de bioinformática
│   └── Aula_Do_FASTQ_a_Taxonomia_16S_Completa.pdf
├── metadata.csv        # Metadados das amostras
├── TUTORIAL.md         # Tutorial completo do pipeline
├── .gitignore
└── README.md           # Este arquivo
```

> Os diretórios `data/`, `databases/`, `scripts/`, `results/`, `logs/` e demais arquivos gerados pela análise **não estão no repositório** (ver [`.gitignore`](.gitignore)). Eles são criados ou linkados localmente conforme descrito no `TUTORIAL.md`.

## Amostras

| ID da Amostra     | Grupo  | Tratamento | Réplica | Descrição                      |
|-------------------|--------|------------|---------|--------------------------------|
| S813_16Sv3v4_01   | liquen | liquen     | 1       | Amostra de liquen              |
| S813_16Sv3v4_02   | liquen | liquen     | 2       | Amostra de liquen              |
| S813_16Sv3v4_03   | musgo  | musgo      | 1       | Amostra de musgo               |
| S813_16Sv3v4_04   | musgo  | musgo      | 2       | Amostra de musgo               |
| S813_16Sv3v4_05   | mel    | mel        | 1       | Amostra de mel                 |
| S813_16Sv3v4_06   | mel    | mel        | 2       | Amostra de mel                 |
| S813_16Sv3v4_07   | mel    | mel        | 3       | Amostra de mel                 |
| S813_16Sv3v4_NN   | —      | Nenhum     | —       | Controle negativo (extração/PCR)|

## Pipeline de Análise

O pipeline completo está descrito em detalhes no arquivo [`TUTORIAL.md`](TUTORIAL.md) (ou [`TUTORIAL.html`](TUTORIAL.html) para visualização renderizada) e inclui:

1. **Controle de qualidade** — FastQC + MultiQC
2. **Remoção de primers** — Cutadapt
3. **Denoising e inferência de ASVs** — DADA2 (modelo de erro *quality-binned*: bins 9, 23, 38)
4. **Atribuição taxonômica** — DADA2 + banco SILVA v138.1
5. **Diversidade alfa** — Índices de riqueza e diversidade (Shannon, Simpson, Observed, Chao1)
6. **Diversidade beta** — Dissimilaridade de Bray-Curtis + PCoA + heatmap + PERMANOVA

## Pré-requisitos

- **R** ≥ 4.3 com os pacotes: `dada2`, `phyloseq`, `vegan`, `ggplot2`, `pheatmap`, `reshape2`, `gridExtra`
- **Cutadapt** 5.2 (via `module load Bio/Cutadapt/5.2`)
- **FastQC** 0.12.1 (via `module load Bio/FastQC/0.12.1`)
- **MultiQC** 1.33 (via `module load Bio/MultiQC/1.33`)
- Acesso ao cluster SLURM para submissão dos scripts bash
- Acesso ao diretório compartilhado do curso no servidor

## Como Começar

1. Clone o repositório e crie os links simbólicos conforme descrito acima
2. Verifique o arquivo `metadata.csv` com as informações das amostras
3. Siga o tutorial passo a passo em `TUTORIAL.md`

## Referências

- Callahan, B.J. et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583.
- Quast, C. et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*, 41, D590–D596.
- McMurdie, P.J. & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE*, 8(4), e61217.

## Licença

Material didático para uso interno no curso de Iniciação Científica em Biotecnologia — LZT0693.
