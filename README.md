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
│   └── slides/
├── parte_02/               # Módulo 2 — Pipeline com novos tratamentos
│   ├── TUTORIAL.md
│   └── metadata.csv
└── parte_03/               # Módulo 3 — Integração e MaAsLin3
    └── TUTORIAL.md
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
* **Tutorial Completo (Integração)**: [`parte_03/TUTORIAL.md`](parte_03/TUTORIAL.md)
* **Tutorial Simplificado (Slice)**: [`parte_03/TUTORIAL_mel_liquen.md`](parte_03/TUTORIAL_mel_liquen.md) - *Ideal para começar! Aborda uma comparação simples de Mel vs Líquen usando apenas a Parte 01.*

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

**2. Compreender o carregamento de módulos no servidor:**

Para que os programas funcionem dentro do servidor, eles já estão instalados, mas precisam ser "chamados" (ou ativados) na sua sessão antes de serem usados. O comando utilizado para ativar um programa é o `module load`, conforme o exemplo abaixo:

```bash
module load Bio/FastQC/0.12.1
module load Bio/MultiQC/1.33
module load Bio/Cutadapt/5.2
```

> **💡 Aviso Importante:** Você **não precisa** executar esses comandos manualmente no seu terminal agora. Esse carregamento será realizado automaticamente dentro dos scripts do SLURM (arquivos `.sh`) que você criará nas próximas etapas. O exemplo acima serve apenas para que você entenda como a ativação de programas funciona na arquitetura do servidor.

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
