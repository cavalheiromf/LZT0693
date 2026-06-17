# LZT0693 — Análise de Amplicon 16S rRNA (V3-V4) | Iniciação Científica em Biotecnologia

## Descrição

Este repositório contém o pipeline de análise de sequenciamento de amplicon **16S rRNA (regiões V3-V4)** desenvolvido na disciplina LZT0693 - Iniciação Científica em Biotecnologia.

O sequenciamento foi realizado em plataforma **Illumina** (paired-end), resultando em 17 bibliotecas: 16 amostras biológicas (Grupo A e Grupo B) e 1 controle negativo.

> **Nota:** Os dados brutos (FASTQ) **não estão incluídos** neste repositório. Eles são acessados por meio de um link simbólico (`data/`) apontando para o diretório compartilhado do servidor.

## Como Clonar / Baixar este Repositório

Para clonar este repositório no terminal do servidor, execute:

```bash
git clone https://github.com/cavalheiromf/LZT0693.git
cd LZT0693
```

Após clonar, recriar os links simbólicos necessários na pasta `parte_02`:

```bash
ln -s /home/cursos/LCoutinho202604/data_shared/aula_03/data     data
ln -s /home/cursos/LCoutinho202604/data_shared/aula_03/databases databases
```

## Estrutura do Repositório (Parte 02)

Arquivos contidos nesta parte:

```
LZT0693/parte_02/
├── metadata.csv        # Metadados das amostras
├── TUTORIAL.md         # Tutorial completo do pipeline
└── README.md           # Este arquivo
```

> Os diretórios `data/`, `databases/`, `scripts/`, `results/`, `logs/` e demais arquivos gerados pela análise **não estão no repositório** (ver [`.gitignore`](.gitignore) na raiz). Eles são criados ou linkados localmente conforme descrito no `TUTORIAL.md`.

## Amostras (Aula 03)

| ID da Amostra     | Grupo    | Tratamento    | Réplica | Descrição                      |
|-------------------|----------|---------------|---------|--------------------------------|
| S813_A1           | Grupo_A  | Tratamento_A  | 1       | Amostra A1                     |
| S813_A2           | Grupo_A  | Tratamento_A  | 2       | Amostra A2                     |
| S813_A3           | Grupo_A  | Tratamento_A  | 3       | Amostra A3                     |
| S813_A4           | Grupo_A  | Tratamento_A  | 4       | Amostra A4                     |
| S813_A5           | Grupo_A  | Tratamento_A  | 5       | Amostra A5                     |
| S813_A6           | Grupo_A  | Tratamento_A  | 6       | Amostra A6                     |
| S813_A7           | Grupo_A  | Tratamento_A  | 7       | Amostra A7                     |
| S813_A8           | Grupo_A  | Tratamento_A  | 8       | Amostra A8                     |
| S813_A9           | Grupo_A  | Tratamento_A  | 9       | Amostra A9                     |
| S813_A10          | Grupo_A  | Tratamento_A  | 10      | Amostra A10                    |
| S813_A11          | Grupo_A  | Tratamento_A  | 11      | Amostra A11                    |
| S813_A12          | Grupo_A  | Tratamento_A  | 12      | Amostra A12                    |
| S813_B1           | Grupo_B  | Tratamento_B  | 1       | Amostra B1                     |
| S813_B2           | Grupo_B  | Tratamento_B  | 2       | Amostra B2                     |
| S813_B3           | Grupo_B  | Tratamento_B  | 3       | Amostra B3                     |
| S813_B4           | Grupo_B  | Tratamento_B  | 4       | Amostra B4                     |
| S813_NN_V3V4      | Controle | Nenhum        | NA      | Controle negativo              |

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
