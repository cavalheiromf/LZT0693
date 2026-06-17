# LZT0693 — Análise de Abundância Diferencial com MaAsLin3 | Parte 03

## Descrição

Este módulo realiza a **integração de dados** de múltiplos sequenciamentos da disciplina LZT0693 (Parte 01 e Parte 02) para avaliar táxons (ASVs) associados a condições biológicas de interesse.

Para isso, utiliza-se a ferramenta **MaAsLin3** (*Microbiome Multivariable Association with Linear Models*), permitindo:
- A mesclagem exata de tabelas de ASVs baseada nas sequências de DNA.
- O ajuste de modelos lineares multivariados.
- O controle estatístico de **efeito de lote (batch effect)** correspondente a cada aula/corrida.

## Estrutura do Repositório (Parte 03)

```
LZT0693/parte_03/
├── README.md           # Este arquivo
├── TUTORIAL.md         # Explicações sobre modelagem e MaAsLin3
└── scripts/
    └── maaslin_analysis.R  # Script principal de fusão de dados e análise R
```

## Como Iniciar

1. Certifique-se de ter concluído o processamento no DADA2 das pastas `parte_01/` e `parte_02/`, pois o script depende da existência de `seqtab_nochim.rds` em ambas.
2. Acesse o RStudio Server.
3. Abra o script [`maaslin_analysis.R`](scripts/maaslin_analysis.R) e execute as etapas para mesclar os dados e rodar o modelo linear.
