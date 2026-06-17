# ==============================================================================
# Análise de Abundância Diferencial com MaAsLin3 (Integração Parte 01 & 02)
# Disciplina LZT0693 - Iniciação Científica em Biotecnologia
# ==============================================================================

# Carregar bibliotecas necessárias
library(dada2)
library(maaslin3)
library(dplyr)
library(phyloseq)

# 1. Configurar caminhos dos dados brutos das aulas anteriores
seqtab1_path <- "caminho/para/parte_01/seqtab_nochim.rds"
seqtab2_path <- "caminho/para/parte_02/seqtab_nochim.rds"

meta1_path <- "caminho/para/parte_01/metadata.csv"
meta2_path <- "caminho/para/parte_02/metadata.csv"

taxa1_path <- "caminho/para/parte_01/taxa.rds"
taxa2_path <- "caminho/para/parte_02/taxa.rds"

# Diretório de saída
output_dir <- "../results/maaslin3"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Verificação se os arquivos das etapas anteriores existem
if (!file.exists(seqtab1_path) || !file.exists(seqtab2_path) ||
    !file.exists(taxa1_path) || !file.exists(taxa2_path)) {
  stop("⚠️ Erro: Tabelas de ASVs ou classificação taxonômica não encontradas em parte_01 ou parte_02. Verifique os caminhos!")
}

cat(">> Carregando tabelas de ASVs e dados taxonômicos...\n")
seqtab1 <- readRDS(seqtab1_path)
seqtab2 <- readRDS(seqtab2_path)
tax1 <- readRDS(taxa1_path)
tax2 <- readRDS(taxa2_path)

cat("Tabela 1 (Parte 01):", nrow(seqtab1), "amostras,", ncol(seqtab1), "ASVs\n")
cat("Tabela 2 (Parte 02):", nrow(seqtab2), "amostras,", ncol(seqtab2), "ASVs\n")

# 2. Mesclar as tabelas de sequências (ASVs idênticas são pareadas automaticamente)
cat("\n>> Mesclando tabelas de ASVs com dada2::mergeSequenceTables...\n")
merged_seqtab <- mergeSequenceTables(seqtab1, seqtab2)
cat("Tabela mesclada final:", nrow(merged_seqtab), "amostras,", ncol(merged_seqtab), "ASVs\n")

# 3. Combinar e limpar a taxonomia das ASVs
cat("\n>> Processando e combinando dados taxonômicos...\n")
combined_tax <- rbind(tax1, tax2)
combined_tax <- combined_tax[!duplicated(rownames(combined_tax)), ]

# 4. Carregar e mesclar os metadados das amostras
cat("\n>> Carregando e processando arquivos de metadados...\n")
meta1 <- read.csv(meta1_path, row.names = 1)
meta2 <- read.csv(meta2_path, row.names = 1)

# Adicionar indicador de lote/batch
meta1$batch <- "parte_01"
meta2$batch <- "parte_02"

# Unir os metadados (garante que as colunas combinem)
merged_metadata <- rbind(meta1, meta2)

# 5. Limpeza de dados (remover controles negativos para análise ecológica/diferencial)
cat("\n>> Filtrando amostras de controle negativo (NN)...\n")
control_samples <- c("S813_16Sv3v4_NN", "S813_NN_V3V4")
merged_metadata <- merged_metadata[!rownames(merged_metadata) %in% control_samples, ]

# Sincronizar tabela de ASVs com os metadados limpos
samples_to_keep <- intersect(rownames(merged_metadata), rownames(merged_seqtab))
merged_seqtab <- merged_seqtab[samples_to_keep, ]
merged_metadata <- merged_metadata[samples_to_keep, ]

# 6. Adicionar a covariável 'reads' (profundidade de sequenciamento/library size)
# Isso permite que o modelo do MaAsLin3 controle estatisticamente variações no volume total de sequências
merged_metadata$reads <- rowSums(merged_seqtab)

# 7. Agrupamento Taxonômico (Taxonomic Glomming) via Phyloseq
cat("\n>> Agrupando abundâncias por nível taxonômico (Família)...\n")
# Criar objeto Phyloseq temporário
ps <- phyloseq(
  otu_table(merged_seqtab, taxa_are_rows = FALSE),
  sample_data(merged_metadata),
  tax_table(combined_tax)
)

# Agrupar táxons no nível de "Family"
# Para agrupar em outros níveis taxonômicos, substitua "Family" por:
# "Kingdom" (Reino), "Phylum" (Filo), "Class" (Classe), "Order" (Ordem), "Genus" (Gênero), ou "Species" (Espécie)
ps_glom <- tax_glom(ps, taxrank = "Family")

# Extrair a matriz de contagens agregada para o MaAsLin3
merged_seqtab_glom <- as(otu_table(ps_glom), "matrix")

# Sincronizar metadados finais contendo a coluna 'reads'
final_metadata <- as(sample_data(ps_glom), "data.frame")

cat("Amostras finais para análise de associação:", nrow(final_metadata), "\n")
cat("Número de táxons agrupados (Família):", ncol(merged_seqtab_glom), "\n")

# 8. Executar o MaAsLin3
# O MaAsLin3 roda modelos lineares generalizados para identificar táxons associados a variáveis de interesse.
# Definimos a fórmula incluindo 'batch' para o efeito de lote e 'reads' para a cobertura do sequenciamento
cat("\n>> Executando MaAsLin3...\n")
maaslin_results <- maaslin3(
  input_data = merged_seqtab_glom,
  input_metadata = final_metadata,
  output = output_dir,
  formula = "~ group + batch + reads",
  normalization = "TSS",    # Total Sum Scaling
  transform = "LOG",        # Transformação logarítmica
  min_abundance = 0.0001,   # Abundância relativa mínima
  min_prevalence = 0.1,     # Prevalência mínima nas amostras
  max_significance = 0.05,  # Limiar de FDR (q-value)
  plot_heatmap = TRUE,
  plot_scatter = TRUE
)

cat("\n>> Análise MaAsLin3 concluída com sucesso!")
cat("\n>> Os resultados foram salvos no diretório:", output_dir, "\n")
