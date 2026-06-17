# ==============================================================================
# Análise de Abundância Diferencial com MaAsLin3
# Disciplina LZT0693 - Iniciação Científica em Biotecnologia
# ==============================================================================

# Configura caminhos das bibliotecas compartilhadas
.libPaths(c("/opt/R/sharedLibs/4.5", .libPaths()))

# Carregar bibliotecas necessárias
library(dada2)
library(maaslin3)
library(dplyr)
library(phyloseq)

setwd('/home/mfcaval/gitlab/aulas/LZT0693/parte_03')

# 1. Configurar caminhos dos arquivos de entrada
seqtab1_path <- "/home/mfcaval/gitlab/aulas/LZT0693/parte_01/results/dada2/seqtab_nochim.rds"
seqtab2_path <- "/home/mfcaval/gitlab/aulas/LZT0693/parte_02/results/dada2/seqtab_nochim.rds"
taxa1_path   <- "/home/mfcaval/gitlab/aulas/LZT0693/parte_01/results/dada2/taxa.rds"
taxa2_path   <- "/home/mfcaval/gitlab/aulas/LZT0693/parte_02/results/dada2/taxa.rds"
meta1_path   <- "/home/mfcaval/gitlab/aulas/LZT0693/parte_01/metadata.csv"
meta2_path   <- "/home/mfcaval/gitlab/aulas/LZT0693/parte_02/metadata.csv"

output_dir   <- "results/maaslin3"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Verificar existência dos arquivos essenciais
if (!all(file.exists(c(seqtab1_path, seqtab2_path, taxa1_path, taxa2_path)))) {
  stop("⚠️ Erro: Arquivos RDS das etapas anteriores não encontrados. Verifique os caminhos!")
}

# 2. Carregar dados brutos
seqtab1 <- readRDS(seqtab1_path)
seqtab2 <- readRDS(seqtab2_path)
tax1    <- readRDS(taxa1_path)
tax2    <- readRDS(taxa2_path)

# 3. Mesclar tabelas de sequências (DADA2)
merged_seqtab <- mergeSequenceTables(seqtab1, seqtab2)
cat("Tabela mesclada:", nrow(merged_seqtab), "amostras,", ncol(merged_seqtab), "ASVs\n")

# 4. Combinar e limpar dados taxonômicos
combined_tax <- rbind(tax1, tax2)
combined_tax <- combined_tax[!duplicated(rownames(combined_tax)), ]

# 5. Carregar e harmonizar metadados
meta1 <- read.csv(meta1_path, row.names = 1)
meta2 <- read.csv(meta2_path, row.names = 1)

# Adicionar covariável de lote (batch effect)
meta1$batch <- "parte_01"
meta2$batch <- "parte_02"
merged_metadata <- rbind(meta1, meta2)

# Remover amostras de controle negativo
control_samples <- c("S813_16Sv3v4_NN", "S813_NN_V3V4")
merged_metadata <- merged_metadata[!rownames(merged_metadata) %in% control_samples, ]

# Sincronizar dados de ASVs e metadados
samples_to_keep <- intersect(rownames(merged_metadata), rownames(merged_seqtab))
merged_seqtab   <- merged_seqtab[samples_to_keep, ]
merged_metadata <- merged_metadata[samples_to_keep, ]

# 5b. [FILTRO DE COMPARAÇÃO OPCIONAL]
# Descomente para comparar apenas grupos/tratamentos específicos de interesse.
# Exemplo (apenas grupo musgo com tratamentos epifita e rupicula):
# merged_metadata <- merged_metadata[
#   merged_metadata$group == "musgo" &
#   merged_metadata$treatment %in% c("epifita", "rupicula"), ]
# samples_to_keep <- intersect(rownames(merged_metadata), rownames(merged_seqtab))
# merged_seqtab   <- merged_seqtab[samples_to_keep, ]
# merged_metadata <- merged_metadata[samples_to_keep, ]

# 6. Calcular a profundidade de sequenciamento
merged_metadata$reads <- rowSums(merged_seqtab)

# 6b. Configurar variáveis categóricas (Fatores)
merged_metadata$treatment <- factor(merged_metadata$treatment)
merged_metadata$group     <- factor(merged_metadata$group)
merged_metadata$batch     <- factor(merged_metadata$batch)

# 6c. Definir níveis de referência (Baselines)
if ("rupicula" %in% levels(merged_metadata$treatment)) {
  merged_metadata$treatment <- relevel(merged_metadata$treatment, ref = "rupicula")
}

# 7. Agrupamento Taxonômico e Construção do Objeto Phyloseq
ps <- phyloseq(
  otu_table(merged_seqtab, taxa_are_rows = FALSE),
  sample_data(merged_metadata),
  tax_table(combined_tax)
)

# Agrupar no nível de Família (substitua por "Genus", "Phylum", etc., se preferir)
ps_glom <- tax_glom(ps, taxrank = "Family")

# Renomear os táxons agrupados (substitui sequências brutas de DNA por nomes legíveis)
tax_mat <- as(tax_table(ps_glom), "matrix")
new_names <- ifelse(
  !is.na(tax_mat[, "Family"]),
  paste(tax_mat[, "Phylum"], tax_mat[, "Family"], sep = "__"),
  ifelse(
    !is.na(tax_mat[, "Order"]),
    paste(tax_mat[, "Phylum"], tax_mat[, "Order"], "unclassified", sep = "__"),
    paste(tax_mat[, "Phylum"], "unclassified", sep = "__")
  )
)
taxa_names(ps_glom) <- make.unique(new_names, sep = "_")

# Extrair tabelas finais sincronizadas para o MaAsLin3
final_seqtab   <- as(otu_table(ps_glom), "matrix")
final_metadata <- as(sample_data(ps_glom), "data.frame")

# 8. Executar o MaAsLin3 com seleção dinâmica de covariáveis
formula_vars <- c("group", "treatment", "batch", "reads")

# Descartar preditores sem variação (evita erros de contrastes com nível único)
vars_to_keep <- sapply(formula_vars, function(v) {
  if (!v %in% colnames(final_metadata)) return(FALSE)
  col <- final_metadata[[v]]
  if (is.factor(col) || is.character(col)) {
    return(length(unique(na.omit(col))) >= 2)
  }
  return(TRUE)
})

vars_kept    <- formula_vars[vars_to_keep]
vars_dropped <- formula_vars[!vars_to_keep]

if (length(vars_dropped) > 0) {
  cat("⚠️ Variáveis removidas da fórmula (apenas 1 nível):", paste(vars_dropped, collapse = ", "), "\n")
}

maaslin_formula <- paste("~", paste(vars_kept, collapse = " + "))
cat("✅ Fórmula estatística:", maaslin_formula, "\n")

maaslin_results <- maaslin3(
  input_data       = final_seqtab,
  input_metadata   = final_metadata,
  output           = output_dir,
  formula          = maaslin_formula,
  normalization    = "TSS",
  transform        = "LOG",
  min_abundance    = 0.0001,
  min_prevalence   = 0.1,
  max_significance = 0.05
)

cat("\n>> Análise concluída! Resultados salvos em:", output_dir, "\n")
