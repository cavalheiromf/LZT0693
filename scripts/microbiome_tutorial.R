# ==============================================================================
# Pipeline Completo de Análise Microbioma 16S rRNA (V3-V4)
# Disciplina LZT0693 - Iniciação Científica em Biotecnologia
#
# Este script foi desenhado para ser executado interativamente no RStudio.
# Ele realiza o processamento das reads trimadas com o Cutadapt,
# executa o pipeline DADA2 (com suporte a qualidade binned do NovaSeq),
# e realiza análises ecológicas básicas de Diversidade Alfa e Beta.
# ==============================================================================

# ==============================================================================
# 1. Carregamento de Pacotes e Instalação (se necessário)
# ==============================================================================

# Se necessário, descomente as linhas abaixo para instalar os pacotes obrigatórios:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("dada2", "phyloseq", "Biostrings"))
# install.packages(c("ggplot2", "vegan", "pheatmap", "reshape2", "gridExtra"))

library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(pheatmap)
library(Biostrings)

# ==============================================================================
# 2. Configuração de Caminhos e Variáveis
# ==============================================================================

# Diretório base do projeto (ajuste se seu diretório ativo for diferente)
base_dir <- "."

# Diretórios de entrada e saída
cut_dir <- file.path(base_dir, "results", "cutadapt")
dada2_dir <- file.path(base_dir, "results", "dada2")

# Criar o diretório de saída para o DADA2 se não existir
if (!dir.exists(dada2_dir)) {
  dir.create(dada2_dir, recursive = TRUE, showWarnings = FALSE)
}

# Localizar os arquivos R1 e R2 trimados produzidos pelo Cutadapt
fnFs <- sort(list.files(cut_dir, pattern = "_R1_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(cut_dir, pattern = "_R2_trimmed.fastq.gz", full.names = TRUE))

# Extrair os nomes das amostras de forma limpa
sample_names <- gsub("_R1_trimmed.fastq.gz", "", basename(fnFs))

cat("--- Verificação inicial ---\n")
cat("Diretório de trabalho:", getwd(), "\n")
cat("Total de amostras encontradas:", length(sample_names), "\n")
print(sample_names)
cat("---------------------------\n\n")

# ==============================================================================
# 3. Inspeção de Perfis de Qualidade
# ==============================================================================
cat(">> Gerando gráficos de perfis de qualidade (Forward/Reverse)...\n")

# Perfil R1 (Forward)
p_fwd <- plotQualityProfile(fnFs, aggregate = TRUE)
ggsave(file.path(dada2_dir, "quality_profile_R1.png"), p_fwd,
       width = 10, height = 6, dpi = 150)

# Perfil R2 (Reverse)
p_rev <- plotQualityProfile(fnRs, aggregate = TRUE)
ggsave(file.path(dada2_dir, "quality_profile_R2.png"), p_rev,
       width = 10, height = 6, dpi = 150)

cat(">> Gráficos de qualidade salvos em:", dada2_dir, "\n\n")

# ==============================================================================
# 4. Filtragem e Trimming (Estratégia MiSeq vs NovaSeq)
# ==============================================================================
cat(">> Iniciando filtragem e trimming das reads...\n")

# Diretório para saída das reads filtradas
filt_dir <- file.path(dada2_dir, "filtered")
if (!dir.exists(filt_dir)) {
  dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
}

# Definir os caminhos dos arquivos de saída filtrados
filtFs <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# ------------------------------------------------------------------------------
# PARÂMETROS - Escolha com base na plataforma (comente/descomente as opções):
# ------------------------------------------------------------------------------
# OPÇÃO A — MiSeq (Phred scores contínuos tradicionais):
#   Geralmente os reads R2 caem muito de qualidade ao final. Ajuste com base nos gráficos.
truncF <- 280
truncR <- 200
truncQ_val <- 2

# OPÇÃO B — NovaSeq (Quality scores binned em degraus - Q37, Q23, Q12, Q2):
#   Remova o comentário das linhas abaixo e comente a Opção A se seus dados forem NovaSeq.
# truncF <- 0
# truncR <- 0
# truncQ_val <- 8
# ------------------------------------------------------------------------------

# Executar o processo de filtragem
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

# Criar tabela de estatísticas de filtragem
filt_stats <- as.data.frame(filt_out)
filt_stats$pct_passed <- round(filt_stats$reads.out / filt_stats$reads.in * 100, 1)

cat("\n>> Estatísticas de filtragem e trimming:\n")
print(filt_stats)
write.csv(filt_stats, file.path(filt_dir, "filter_stats.csv"), row.names = TRUE)
cat("\n")

# ==============================================================================
# 5. Aprendizado do Modelo de Erros (Binned Quality Support)
# ==============================================================================
cat(">> Aprendendo modelos de taxas de erro com DADA2...\n")

# Configurar função para lidar com bins discretos de qualidade (12, 24, 40)
# essencial para dados NovaSeq e robusto para demais plataformas
binned_err_fun <- makeBinnedQualErrfun(nqual = 3, binnedQuals = c(12, 24, 40))

# Aprender erro R1
cat(">> Aprendendo erros do Forward (R1)...\n")
errF <- learnErrors(filtFs, errorEstimationFunction = binned_err_fun, multithread = TRUE)

# Aprender erro R2
cat(">> Aprendendo erros do Reverse (R2)...\n")
errR <- learnErrors(filtRs, errorEstimationFunction = binned_err_fun, multithread = TRUE)

# Alternativa MiSeq padrão (descomente caso queira usar a regressão loess padrão):
# errF <- learnErrors(filtFs, multithread = TRUE)
# errR <- learnErrors(filtRs, multithread = TRUE)

# Gerar e salvar gráficos do modelo de erro
p_errF <- plotErrors(errF, nominalQ = TRUE)
ggsave(file.path(dada2_dir, "error_model_Forward.png"), p_errF, width = 10, height = 8, dpi = 150)

p_errR <- plotErrors(errR, nominalQ = TRUE)
ggsave(file.path(dada2_dir, "error_model_Reverse.png"), p_errR, width = 10, height = 8, dpi = 150)

cat(">> Gráficos dos modelos de erro salvos.\n\n")

# ==============================================================================
# 6. Inferência de Variantes de Sequência de Amplicon (ASVs)
# ==============================================================================
cat(">> Executando inferência de ASVs de alta resolução (Denoising)...\n")

# Inferir sequências reais Forward e Reverse
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

cat("\nResumo de ASVs exclusivas inferidas por amostra (Forward):\n")
print(sapply(dadaFs, function(x) length(x$denoised)))
cat("\n")

# ==============================================================================
# 7. Merge de Reads Paired-End e Remoção de Quimeras
# ==============================================================================
cat(">> Realizando o alinhamento e merge de R1 e R2 pareados...\n")
merged <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Construir tabela de frequências de sequências (seqtab)
seqtab <- makeSequenceTable(merged)
cat("\nDimensões iniciais da tabela de ASVs (Amostras x ASVs):", dim(seqtab), "\n")
cat("Distribuição do tamanho em pares de bases dos amplicons:\n")
print(table(nchar(getSequences(seqtab))))

# Remover chimeras bimeras (falsos amplicons formados na PCR)
cat("\n>> Executando detecção e remoção de quimeras...\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

cat("\nASVs retidas após quimeras:", ncol(seqtab_nochim), "de", ncol(seqtab), "\n")
cat("Porcentagem de reads retidos:", round(sum(seqtab_nochim) / sum(seqtab) * 100, 1), "%\n\n")

# ==============================================================================
# 8. Tracking (Rastreamento de Reads por Etapa)
# ==============================================================================
cat(">> Compilando tracking de reads ao longo do pipeline...\n")

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
cat(">> Iniciando atribuição taxonômica com SILVA v138.2...\n")

# Caminhos do banco (certifique-se de baixar nas pastas corretas antes)
silva_train <- file.path(base_dir, "databases", "silva_nr99_v138.2_train_set.fa.gz")
silva_species <- file.path(base_dir, "databases", "silva_species_assignment_v138.2.fa.gz")

if (!file.exists(silva_train) || !file.exists(silva_species)) {
  cat("⚠️ ALERTA: Banco de dados SILVA não encontrado em 'databases/'.\n")
  cat("Por favor, realize o download dos arquivos em: https://zenodo.org/records/4587955\n")
  cat("Interrompendo a execução automática da taxonomia. Salve os dados intermediários e execute depois.\n")
} else {
  # Atribuição taxonômica até gênero
  cat("Atribuindo classificação taxonômica (até nível de Gênero)...\n")
  taxa <- assignTaxonomy(seqtab_nochim, silva_train, multithread = TRUE)
  
  # Atribuição taxonômica de espécie (exata)
  cat("Realizando atribuição taxonômica exata de Espécie...\n")
  taxa <- addSpecies(taxa, silva_species)
  
  cat("\nVisualização das primeiras classificações:\n")
  taxa_print <- taxa
  rownames(taxa_print) <- NULL  # Ocultar a sequência longa de DNA na tela
  print(head(taxa_print))
  
  # Salvar resultados do DADA2
  saveRDS(seqtab_nochim, file.path(dada2_dir, "seqtab_nochim.rds"))
  saveRDS(taxa, file.path(dada2_dir, "taxa.rds"))
  write.csv(taxa, file.path(dada2_dir, "taxonomy_table.csv"))
  cat(">> Tabelas de taxonomia salvas em:", dada2_dir, "\n\n")
}

# ==============================================================================
# 10. Construção do Objeto Phyloseq
# ==============================================================================
cat(">> Estruturando objeto do phyloseq...\n")

metadata_file <- file.path(base_dir, "metadata.csv")

if (!file.exists(metadata_file)) {
  cat("⚠️ ALERTA: Arquivo 'metadata.csv' não encontrado. Por favor, crie-o antes de prosseguir.\n")
} else if (exists("taxa")) {
  # Carregar os metadados do experimento
  metadata <- read.csv(metadata_file, row.names = 1)
  
  # Criar objeto unificado phyloseq
  ps <- phyloseq(
      otu_table(seqtab_nochim, taxa_are_rows = FALSE),
      sample_data(metadata),
      tax_table(taxa)
  )
  
  # Renomear as sequências genômicas longas dos nomes das taxa para ASV1, ASV2...
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  
  cat("\nObjeto phyloseq criado com sucesso:\n")
  print(ps)
  saveRDS(ps, file.path(dada2_dir, "phyloseq_object.rds"))
  cat("\n")
}

# ==============================================================================
# 11. Análise Ecológica - Diversidade Alfa
# ==============================================================================
if (exists("ps")) {
  cat(">> Calculando métricas de Diversidade Alfa...\n")
  
  # 11.1 Filtragens ecológicas básicas (remover controle negativo e ASVs com abundância 0)
  ps_bio <- prune_samples(sample_names(ps) != "S813_16Sv3v4_NN", ps)
  ps_bio <- prune_taxa(taxa_sums(ps_bio) > 0, ps_bio)
  
  # Estimar métricas de riqueza e diversidade
  alpha_div <- estimate_richness(ps_bio, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
  alpha_div$sample_id <- rownames(alpha_div)
  
  cat("\nTabela de Diversidade Alfa calculada:\n")
  print(alpha_div)
  write.csv(alpha_div, file.path(dada2_dir, "alpha_diversity.csv"), row.names = FALSE)
  
  # 11.2 Plotar e salvar gráfico de diversidade alfa
  p_obs <- plot_richness(ps_bio, measures = c("Observed", "Shannon", "Simpson")) +
      geom_boxplot(alpha = 0.3) +
      theme_bw() +
      theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          strip.text = element_text(size = 12, face = "bold")
      ) +
      labs(title = "Diversidade Alfa — 16S V3-V4")
  
  ggsave(file.path(dada2_dir, "alpha_diversity_plot.png"), p_obs, width = 12, height = 6, dpi = 150)
  cat(">> Gráfico de diversidade alfa salvo em 'alpha_diversity_plot.png'.\n\n")
  
  # ----------------------------------------------------------------------------
  # Se quiser plotar por grupo experimental (após editar metadata.csv):
  # ----------------------------------------------------------------------------
  # p_group <- plot_richness(ps_bio, x = "group", measures = c("Observed", "Shannon", "Simpson"), color = "group") +
  #     geom_boxplot(alpha = 0.3) +
  #     theme_bw() +
  #     labs(title = "Diversidade Alfa por Grupo")
  # ggsave(file.path(dada2_dir, "alpha_diversity_by_group.png"), p_group, width = 12, height = 6, dpi = 150)
}

# ==============================================================================
# 12. Análise Ecológica - Diversidade Beta (Bray-Curtis)
# ==============================================================================
if (exists("ps_bio")) {
  cat(">> Calculando métricas de Diversidade Beta (Bray-Curtis)...\n")
  
  # Normalizar abundância (converter abundância absoluta em abundância relativa)
  ps_rel <- transform_sample_counts(ps_bio, function(x) x / sum(x))
  
  # 12.1 Calcular matriz de dissimilaridade de Bray-Curtis
  bray_dist <- phyloseq::distance(ps_rel, method = "bray")
  bray_mat <- as.matrix(bray_dist)
  
  cat("\nMatriz de dissimilaridade Bray-Curtis:\n")
  print(bray_mat)
  write.csv(bray_mat, file.path(dada2_dir, "bray_curtis_matrix.csv"))
  
  # 12.2 Ordenação e PCoA
  pcoa_res <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)
  
  p_pcoa <- plot_ordination(ps_rel, pcoa_res, type = "samples") +
      geom_point(size = 4, alpha = 0.8) +
      theme_bw() +
      theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12)
      ) +
      labs(title = "PCoA — Bray-Curtis — 16S V3-V4") +
      stat_ellipse(level = 0.95, linetype = 2)
  
  ggsave(file.path(dada2_dir, "pcoa_bray_curtis.png"), p_pcoa, width = 8, height = 6, dpi = 150)
  
  # 12.3 Heatmap da matriz de Bray-Curtis
  pheatmap(bray_mat,
           clustering_distance_rows = bray_dist,
           clustering_distance_cols = bray_dist,
           color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
           main = "Bray-Curtis — Heatmap",
           fontsize = 10,
           filename = file.path(dada2_dir, "bray_curtis_heatmap.png"),
           width = 8, height = 7)
  
  # 12.4 PERMANOVA (teste estatístico)
  # (Descomente para executar quando preencher o metadata.csv com seus grupos biológicos reais)
  # meta_df <- data.frame(sample_data(ps_bio))
  # set.seed(42)
  # permanova_res <- adonis2(bray_dist ~ group, data = meta_df, permutations = 999)
  # cat("\n>> Resultado PERMANOVA:\n")
  # print(permanova_res)
  # write.csv(as.data.frame(permanova_res), file.path(dada2_dir, "permanova_bray_curtis.csv"))
  
  cat("\n>> Análises ecológicas concluídas! Todos os plots e tabelas foram salvos no diretório dada2.\n")
}
