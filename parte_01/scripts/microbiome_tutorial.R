# ==============================================================================
# Pipeline Completo de Análise Microbioma 16S rRNA (V3-V4)
# Disciplina LZT0693 - Iniciação Científica em Biotecnologia
#
# Este script realiza o processamento das reads trimadas com o Cutadapt,
# executa o pipeline DADA2 e realiza análises ecológicas.
# ==============================================================================

library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(pheatmap)
library(Biostrings)

# Função auxiliar para log de tempo
log_time <- function(step_name, start_time) {
  end_time <- Sys.time()
  diff <- difftime(end_time, start_time, units = "secs")
  cat(sprintf("\n>> [TEMPO] %s concluído em %.2f segundos (%.2f minutos).\n\n", step_name, diff, diff/60))
  return(end_time)
}

start_all <- Sys.time()

# ==============================================================================
# 2. Configuração de Caminhos e Variáveis
# ==============================================================================
base_dir <- "."
cut_dir <- file.path(base_dir, "results", "cutadapt")
dada2_dir <- file.path(base_dir, "results", "dada2")

if (!dir.exists(dada2_dir)) {
  dir.create(dada2_dir, recursive = TRUE, showWarnings = FALSE)
}

fnFs <- sort(list.files(cut_dir, pattern = "_R1_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(cut_dir, pattern = "_R2_trimmed.fastq.gz", full.names = TRUE))
sample_names <- gsub("_R1_trimmed.fastq.gz", "", basename(fnFs))

cat("--- Verificação inicial ---\n")
cat("Diretório de trabalho:", getwd(), "\n")
cat("Total de amostras encontradas:", length(sample_names), "\n")
print(sample_names)
cat("---------------------------\n\n")

# ==============================================================================
# 3. Inspeção de Perfis de Qualidade
# ==============================================================================
t_step <- Sys.time()
cat(">> Gerando gráficos de perfis de qualidade (Forward/Reverse)...\n")

p_fwd <- plotQualityProfile(fnFs, aggregate = TRUE)
ggsave(file.path(dada2_dir, "quality_profile_R1.png"), p_fwd, width = 10, height = 6, dpi = 150)

p_rev <- plotQualityProfile(fnRs, aggregate = TRUE)
ggsave(file.path(dada2_dir, "quality_profile_R2.png"), p_rev, width = 10, height = 6, dpi = 150)

t_step <- log_time("Inspeção de Perfis de Qualidade", t_step)

# ==============================================================================
# 4. Filtragem e Trimming
# ==============================================================================
cat(">> Iniciando filtragem e trimming das reads...\n")

filt_dir <- file.path(dada2_dir, "filtered")
if (!dir.exists(filt_dir)) {
  dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
}

filtFs <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# Parâmetros MiSeq (dados de qualidade contínuos)
truncF <- 280
truncR <- 200
truncQ_val <- 2

filt_out <- filterAndTrim(
    fnFs, filtFs,
    fnRs, filtRs,
    truncLen = c(truncF, truncR),
    maxN = 0,
    maxEE = c(2, 2),
    truncQ = truncQ_val,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
)

filt_stats <- as.data.frame(filt_out)
filt_stats$pct_passed <- round(filt_stats$reads.out / filt_stats$reads.in * 100, 1)

cat("\n>> Estatísticas de filtragem e trimming:\n")
print(filt_stats)
write.csv(filt_stats, file.path(filt_dir, "filter_stats.csv"), row.names = TRUE)

t_step <- log_time("Filtragem e Trimming", t_step)

# ==============================================================================
# 5. Aprendizado do Modelo de Erros
# ==============================================================================
cat(">> Aprendendo modelos de taxas de erro com DADA2...\n")

cat(">> Aprendendo erros do Forward (R1)...\n")
errF <- learnErrors(filtFs, multithread = TRUE)

cat(">> Aprendendo erros do Reverse (R2)...\n")
errR <- learnErrors(filtRs, multithread = TRUE)

p_errF <- plotErrors(errF, nominalQ = TRUE)
ggsave(file.path(dada2_dir, "error_model_Forward.png"), p_errF, width = 10, height = 8, dpi = 150)

p_errR <- plotErrors(errR, nominalQ = TRUE)
ggsave(file.path(dada2_dir, "error_model_Reverse.png"), p_errR, width = 10, height = 8, dpi = 150)

t_step <- log_time("Aprendizado do Modelo de Erros", t_step)

# ==============================================================================
# 6. Inferência de Variantes de Sequência de Amplicon (ASVs) - Denoising
# ==============================================================================
cat(">> Executando inferência de ASVs de alta resolução (Denoising)...\n")

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

cat("\nResumo de ASVs exclusivas inferidas por amostra (Forward):\n")
print(sapply(dadaFs, function(x) length(x$denoised)))

t_step <- log_time("Inferência de ASVs (Denoising)", t_step)

# ==============================================================================
# 7. Merge de Reads Paired-End e Remoção de Quimeras
# ==============================================================================
cat(">> Realizando o alinhamento e merge de R1 e R2 pareados...\n")
merged <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab <- makeSequenceTable(merged)
cat("\nDimensões iniciais da tabela de ASVs (Amostras x ASVs):", dim(seqtab), "\n")
cat("Distribuição do tamanho em pb dos amplicons:\n")
print(table(nchar(getSequences(seqtab))))

cat("\n>> Executando detecção e remoção de quimeras...\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

cat("\nASVs retidas após quimeras:", ncol(seqtab_nochim), "de", ncol(seqtab), "\n")
cat("Porcentagem de reads retidos:", round(sum(seqtab_nochim) / sum(seqtab) * 100, 1), "%\n\n")

t_step <- log_time("Merge e Remoção de Quimeras", t_step)

# ==============================================================================
# 8. Tracking (Rastreamento de Reads por Etapa)
# ==============================================================================
getN <- function(x) sum(getUniques(x))

track <- data.frame(
    input      = filt_out[, 1],
    filtered   = filt_out[, 2],
    denoisedF  = sapply(dadaFs, getN),
    denoisedR  = sapply(dadaRs, getN),
    merged     = sapply(merged, getN),
    nonchim    = rowSums(seqtab_nochim)
)
rownames(track) <- sample_names

cat("\n--- Tabela de Rastreamento (Reads) ---\n")
print(track)
write.csv(track, file.path(dada2_dir, "read_tracking.csv"))
cat("----------------------------------------\n\n")

# ==============================================================================
# 9. Atribuição Taxonômica (Banco SILVA)
# ==============================================================================
cat(">> Iniciando atribuição taxonômica com SILVA v138.1...\n")

silva_train <- file.path(base_dir, "databases", "silva_nr99_v138.1_train_set.fa.gz")
silva_species <- file.path(base_dir, "databases", "silva_species_assignment_v138.1.fa.gz")

if (!file.exists(silva_train) || !file.exists(silva_species)) {
  cat("⚠️ ALERTA: Banco de dados SILVA não encontrado em 'databases/'.\n")
  cat("Por favor, verifique se o link simbólico para a pasta 'databases' compartilhada no servidor foi criado corretamente.\n")
} else {
  cat("Atribuindo classificação taxonômica (até nível de Gênero)...\n")
  taxa <- assignTaxonomy(seqtab_nochim, silva_train, multithread = TRUE)
  
  cat("Realizando atribuição taxonômica exata de Espécie...\n")
  taxa <- addSpecies(taxa, silva_species)
  
  cat("\nVisualização das primeiras classificações:\n")
  taxa_print <- taxa
  rownames(taxa_print) <- NULL
  print(head(taxa_print))
  
  saveRDS(seqtab_nochim, file.path(dada2_dir, "seqtab_nochim.rds"))
  saveRDS(taxa, file.path(dada2_dir, "taxa.rds"))
  write.csv(taxa, file.path(dada2_dir, "taxonomy_table.csv"))
  cat(">> Tabelas de taxonomia salvas em:", dada2_dir, "\n\n")
}

t_step <- log_time("Atribuição Taxonômica (SILVA)", t_step)

# ==============================================================================
# 10. Construção do Objeto Phyloseq
# ==============================================================================
cat(">> Estruturando objeto do phyloseq...\n")

metadata_file <- file.path(base_dir, "metadata.csv")

if (!file.exists(metadata_file)) {
  cat("⚠️ ALERTA: Arquivo 'metadata.csv' não encontrado.\n")
} else if (exists("taxa")) {
  metadata <- read.csv(metadata_file, row.names = 1)
  
  ps <- phyloseq(
      otu_table(seqtab_nochim, taxa_are_rows = FALSE),
      sample_data(metadata),
      tax_table(taxa)
  )
  
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  
  cat("\nObjeto phyloseq criado com sucesso:\n")
  print(ps)
  saveRDS(ps, file.path(dada2_dir, "phyloseq_object.rds"))
}

t_step <- log_time("Construção do Objeto Phyloseq", t_step)

# ==============================================================================
# 11. Análise Ecológica - Diversidade Alfa
# ==============================================================================
if (exists("ps")) {
  cat(">> Calculando métricas de Diversidade Alfa...\n")
  
  ps_bio <- prune_samples(sample_names(ps) != "S813_16Sv3v4_NN", ps)
  ps_bio <- prune_taxa(taxa_sums(ps_bio) > 0, ps_bio)
  
  alpha_div <- estimate_richness(ps_bio, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
  alpha_div$sample_id <- rownames(alpha_div)
  
  cat("\nTabela de Diversidade Alfa calculada:\n")
  print(alpha_div)
  write.csv(alpha_div, file.path(dada2_dir, "alpha_diversity.csv"), row.names = FALSE)
  
  p_obs <- plot_richness(ps_bio, x = "group", measures = c("Observed", "Shannon", "Simpson"), color = "group") +
      geom_boxplot(alpha = 0.3) +
      theme_bw() +
      theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          strip.text = element_text(size = 12, face = "bold")
      ) +
      labs(title = "Diversidade Alfa por Grupo — 16S V3-V4")
  
  ggsave(file.path(dada2_dir, "alpha_diversity_plot.png"), p_obs, width = 12, height = 6, dpi = 150)
  cat(">> Gráfico de diversidade alfa salvo.\n\n")
}

t_step <- log_time("Diversidade Alfa", t_step)

# ==============================================================================
# 12. Análise Ecológica - Diversidade Beta (Bray-Curtis) e PERMANOVA
# ==============================================================================
if (exists("ps_bio")) {
  cat(">> Calculando métricas de Diversidade Beta (Bray-Curtis)...\n")
  
  ps_rel <- transform_sample_counts(ps_bio, function(x) x / sum(x))
  
  bray_dist <- phyloseq::distance(ps_rel, method = "bray")
  bray_mat <- as.matrix(bray_dist)
  
  cat("\nMatriz de dissimilaridade Bray-Curtis:\n")
  print(bray_mat)
  write.csv(bray_mat, file.path(dada2_dir, "bray_curtis_matrix.csv"))
  
  pcoa_res <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)
  
  p_pcoa <- plot_ordination(ps_rel, pcoa_res, type = "samples", color = "group") +
      geom_point(size = 4, alpha = 0.8) +
      theme_bw() +
      theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12)
      ) +
      labs(title = "PCoA — Bray-Curtis — 16S V3-V4") +
      stat_ellipse(level = 0.95, linetype = 2)
  
  ggsave(file.path(dada2_dir, "pcoa_bray_curtis.png"), p_pcoa, width = 8, height = 6, dpi = 150)
  
  pheatmap(bray_mat,
           clustering_distance_rows = bray_dist,
           clustering_distance_cols = bray_dist,
           color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
           main = "Bray-Curtis — Heatmap",
           fontsize = 10,
           filename = file.path(dada2_dir, "bray_curtis_heatmap.png"),
           width = 8, height = 7)
  
  cat(">> Executando PERMANOVA...\n")
  meta_df <- data.frame(sample_data(ps_bio))
  set.seed(42)
  permanova_res <- adonis2(bray_dist ~ group, data = meta_df, permutations = 999)
  cat("\n>> Resultado PERMANOVA:\n")
  print(permanova_res)
  write.csv(as.data.frame(permanova_res), file.path(dada2_dir, "permanova_bray_curtis.csv"))
}

t_step <- log_time("Diversidade Beta e PERMANOVA", t_step)

# ==============================================================================
# 13. Composição Taxonômica (Abundância Relativa)
# ==============================================================================
if (exists("ps_bio")) {
  cat(">> Gerando gráficos de Composição Taxonômica por Filo...\n")
  
  # Agrupar/aglomerar por nível de Filo
  ps_phylum <- tax_glom(ps_bio, "Phylum")
  
  # Converter abundância absoluta para abundância relativa
  ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
  
  # Plotar composição taxonômica (gráfico de barras)
  p_bar <- plot_bar(ps_phylum_rel, fill = "Phylum") +
      geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
      theme_bw() +
      theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(size = 14, face = "bold")
      ) +
      labs(title = "Composição Taxonômica por Filo (Abundância Relativa)",
           x = "Amostra",
           y = "Abundância Relativa")
  
  ggsave(file.path(dada2_dir, "taxonomic_composition_phylum.png"), p_bar, width = 10, height = 6, dpi = 150)
  cat(">> Gráfico de abundância relativa salvo em 'taxonomic_composition_phylum.png'.\n\n")
}

t_step <- log_time("Composição Taxonômica", t_step)

total_diff <- difftime(Sys.time(), start_all, units = "secs")
cat(sprintf(">> [TEMPO] PIPELINE COMPLETO concluído em %.2f segundos (%.2f minutos).\n\n", total_diff, total_diff/60))
