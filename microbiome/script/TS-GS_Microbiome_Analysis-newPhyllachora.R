# Host-Specific Fungal Communities and Early Infection Dynamics of
# Tar Spot on Goosegrass and Corn
#
# Authors: Larissa C. Ferreira, Vitor A. S. Moura, Fernanda R. Silva,
#          Katia V. Xavier
# Affiliation: Everglades Research and Education Center,
#              University of Florida/IFAS, Belle Glade FL 33430, USA
# Contact: kvianaxavier@ufl.edu
#
# Description:
#   Full microbiome analysis pipeline for ITS amplicon sequencing data
#   from tar spot lesions on corn (Zea mays) and goosegrass (Eleusine indica).
#   Covers: DADA2 denoising, phyloseq object construction, normalisation,
#   alpha and beta diversity, community dispersion, Venn/core microbiome
#   analysis, phylogenetic tree visualisation, differential abundance
#   analysis, indicator species analysis, and Phyllachora identity and
#   pathobiome-centered association analyses.
#
# Input:
#   - Paired-end Illumina reads in SRA PRJNA1511899
#   - UNITE ITS database (sh_general_release_dynamic_19.02.2025.fasta)
#
#
# R version: >= 4.3.0



# 1. Load packages ####

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  # DADA2 pipeline
  dada2, ShortRead, Biostrings,
  # Phyloseq and community ecology
  phyloseq, vegan, MicEco, indicspecies, pwalign,
  # Phylogenetics
  DECIPHER, phangorn, ape, pegas, msa,ggtree, parallel,
  # Differential abundance and normalisation
  edgeR,
  # Statistical modelling
  lme4, lmerTest, glmmTMB, emmeans, multcomp, ordinal, broom,
  # Visualisation
  tidyverse, ggpubr, ggthemes, ggrepel, grid, DHARMa, ggResidpanel,
  patchwork, ComplexUpset, ggupset, scales)

set.seed(123)

# 2. Publication theme ####

theme_Publication <- function(base_size = 12, base_family = "Arial") {
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title       = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      panel.border     = element_blank(),
      axis.title       = element_text(face = "bold", size = rel(1)),
      axis.title.y     = element_text(angle = 90, vjust = 2),
      axis.title.x     = element_text(vjust = -0.2),
      axis.line        = element_line(colour = "black"),
      axis.ticks       = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key       = element_blank(),
      legend.position  = "right",
      legend.direction = "vertical",
      legend.key.size  = unit(0.2, "cm"),
      legend.margin    = unit(0.5, "cm"),
      legend.justification = "center",
      plot.margin      = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text       = element_text(face = "bold")
    )
}

host_colours <- c(
  "Corn PL"       = "#648FFF",
  "Corn SL"       = "#785EF0",
  "Goosegrass PL" = "#DC267F",
  "Goosegrass SL" = "#FE6100"
)

host_colours_host <- c(
  "Corn"       = "#648FFF",
  "Goosegrass" = "#FE6100"
)

# 3. DADA2: quality filtering, denoising, and ASV inference ####

path <- "../Reads"

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz",
                        recursive = TRUE, full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz",
                        recursive = TRUE, full.names = TRUE))

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)

sample.names <- sapply(substr(basename(fnFs), 1, 11), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  minLen = 50, maxN = 0, maxEE = c(2, 2),
  truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = TRUE
)

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab        <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)

cat("Fraction of reads retained after chimera removal:",
    round(sum(seqtab.nochim) / sum(seqtab), 3), "\n")

getN  <- function(x) sum(getUniques(x))
track <- cbind(
  out,
  sapply(dadaFs,  getN),
  sapply(dadaRs,  getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, "reads_track.csv")


# 4. Phyloseq object construction ####

## 4a. Taxonomy assignment ####
unite.ref <- "microbiome/data/UNITE_19.02.2025_plus_new_Phyllachora.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref,
                       multithread = TRUE, tryRC = TRUE, verbose = TRUE)

## 4b. Sample metadata ####
# Sample name structure: characters 1-6 = sample ID, 8-9 = lesion type,
# 11 = replicate number
samples.out <- rownames(seqtab.nochim)

samdf <- data.frame(
  SampleID = samples.out,
  samples  = samples.out,
  Sample   = substr(samples.out, 1, 6),
  rep      = substr(samples.out, 11, 11),
  Lesion   = substr(samples.out, 8, 9),
  stringsAsFactors = FALSE
)

rownames(samdf) <- samdf$SampleID

samdf
samdf$rep <- ifelse(samdf$rep == "S", "1", samdf$rep)
samdf$Host <- ifelse(samdf$Sample == "VM1247", "Corn", "Goosegrass")
samdf$Location <- ifelse(samdf$Sample %in% c("VM1245", "VM1247"), "2", "1")
samdf$Group <- paste0(samdf$Host, " ", samdf$Location)
samdf$Group_lesion <- paste0(samdf$Host, " ", samdf$Lesion)
samdf$Sample_lesion <- paste0(samdf$Sample, " ", samdf$Lesion)

rownames(samdf) <- samdf$SampleID


## 4c. Reference sequences ####
ASVs.nochim        <- DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) <- paste0("ASV", seq_len(ncol(seqtab.nochim)))

## 4d. Phylogenetic tree (neighbour-joining + GTR ML refinement) ####
alignment   <- AlignSeqs(ASVs.nochim, anchor = NA, processors = 30)
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm          <- dist.ml(phang.align)
treeNJ      <- NJ(dm)
fit         <- pml(treeNJ, data = phang.align)
fitGTR      <- update(fit, k = 4, inv = 0.2)

# Maximum-likelihood optimization under GTR + Gamma + invariant sites
fitGTR <- optim.pml(
  fit,
  model = "GTR",
  optInv = TRUE,
  optGamma = TRUE,
  optNni = TRUE,
  rearrangement = "stochastic",
  control = pml.control(trace = 1)
)

# Bootstrap support
set.seed(123)

# Bootstrap ML tree with visible progress ####
n_boot  <- 500
batch_n <- 50
n_cores <- 4

n_batches <- ceiling(n_boot / batch_n)

bs_list <- vector("list", n_batches)

start_time <- Sys.time()

for (i in seq_len(n_batches)) {
  
  cat(
    sprintf(
      "[%s] Starting batch %d/%d (%d bootstrap replicates)\n",
      format(Sys.time(), "%H:%M:%S"),
      i,
      n_batches,
      batch_n
    )
  )
  
  batch_start <- Sys.time()
  
  bs_list[[i]] <- bootstrap.pml(
    fitGTR,
    bs = batch_n,
    optNni = TRUE,
    multicore = TRUE,
    mc.cores = n_cores
  )
  
  batch_time <- difftime(
    Sys.time(),
    batch_start,
    units = "mins"
  )
  
  completed <- i * batch_n
  
  cat(
    sprintf(
      "Completed %d/%d bootstraps (%.1f%%) | batch time: %.1f min\n\n",
      completed,
      n_boot,
      100 * completed / n_boot,
      as.numeric(batch_time)
    )
  )
  
  # Save after every batch
  saveRDS(
    bs_list,
    "bootstrap_batches_partial.rds"
  )
}

bs <- do.call(c, bs_list)

class(bs) <- "multiPhylo"

class(bs)

length(bs)

# Add bootstrap support to ML tree
tree_boot <- plotBS(
  fitGTR$tree,
  bs,
  type = "none"
)

all_tips_match <- all(
  sapply(
    bs,
    function(x) setequal(x$tip.label, fitGTR$tree$tip.label)
  )
)

all_tips_match

## 4e. Format taxonomy table ####
tmp.seqtab           <- seqtab.nochim
colnames(tmp.seqtab) <- names(ASVs.nochim)
rownames(taxa)       <- names(ASVs.nochim)

prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
for (prefix in prefixes) taxa <- gsub(prefix, "", taxa)
taxa[is.na(taxa)] <- "Uncharacterized"

## 4f. Assemble phyloseq object ####
ps <- phyloseq(
  otu_table(tmp.seqtab, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa),
  refseq(ASVs.nochim),
  phy_tree(tree_boot)
)
ps

stopifnot(identical(sample_names(ps), rownames(sample_data(ps))))
stopifnot(all(sample_names(ps) == samdf$SampleID))
print(data.frame(sample_data(ps)) %>% count(Host, Group_lesion))

# 5. Pre-processing and filtering ####

# Remove singletons and taxa without phylum assignment
ps <- prune_taxa(taxa_sums(ps) > 1, ps)
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "Uncharacterized"))

together  <- psmelt(ps)
together1 <- together %>%
  filter(OTU != "NA1")

# Retain only ASVs with >= 0.1% relative abundance in at least one sample
# (Nikodemova et al. 2023, Front Cell Infect Microbiol)
together2 <- together1 %>%
  group_by(Sample) %>%
  mutate(
    total_reads = sum(Abundance),
    percent     = Abundance / total_reads * 100
  ) %>%
  ungroup() %>%
  filter(percent >= 0.1)

ps1 <- prune_taxa(taxa_names(ps) %in% together2$OTU, ps)

stopifnot(identical(sample_names(ps1), rownames(sample_data(ps1))))
print(data.frame(sample_data(ps1)) %>% count(Host, Group_lesion))

saveRDS(ps1, "ps1.rds")


# 6. Normalisation with edgeR (relative log expression) ####

edgeRnorm <- function(physeq, ...) {
  require("edgeR")
  require("phyloseq")
  
  if (!taxa_are_rows(physeq)) physeq <- t(physeq)
  
  x <- as(otu_table(physeq), "matrix")
  storage.mode(x) <- "numeric"
  
  y <- edgeR::DGEList(counts = x, remove.zeros = TRUE)
  z <- edgeR::calcNormFactors(y, ...)
  
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("edgeR::calcNormFactors returned non-finite norm factors.")
  }
  
  return(z)
}

tax_table_original   <- tax_table(ps1)
sample_data_original <- sample_data(ps1)
phy_tree_original    <- phy_tree(ps1)

write.csv(tax_table_original, "taxtable.csv", row.names = TRUE, quote = FALSE)

ps1_norm       <- edgeRnorm(ps1)
ps_norm_counts <- edgeR::equalizeLibSizes(ps1_norm)

norm_mat <- ps_norm_counts$pseudo.counts
storage.mode(norm_mat) <- "numeric"

# Integrity checks: confirm normalised matrix matches ps1 exactly before
# rebuilding the phyloseq object, since edgeR can silently drop all-zero
# samples/taxa during DGEList construction.
stopifnot(all(sample_names(ps1) %in% colnames(norm_mat)))
stopifnot(all(taxa_names(ps1) %in% rownames(norm_mat)))

norm_mat <- norm_mat[taxa_names(ps1), sample_names(ps1), drop = FALSE]

ps_norm <- phyloseq(
  otu_table(norm_mat, taxa_are_rows = TRUE),
  tax_table_original,
  sample_data_original,
  phy_tree_original
)

stopifnot(identical(sample_names(ps_norm), rownames(sample_data(ps_norm))))
print(data.frame(sample_data(ps_norm)) %>% count(Host, Group_lesion))

saveRDS(ps_norm, "ps_norm.rds")


# 7. Alpha diversity ####

otu_matrix <- as.matrix(otu_table(ps_norm))

shannon_index <- apply(
  otu_matrix,
  2,
  function(x) vegan::diversity(x, index = "shannon")
)

simpson_index <- apply(
  otu_matrix,
  2,
  function(x) vegan::diversity(x, index = "simpson")
)
observed_richness <- apply(otu_matrix, 2, function(x) sum(x > 0))

richness_diversity_manual <- data.frame(
  SampleID = colnames(otu_matrix),
  Observed = observed_richness,
  Shannon  = shannon_index,
  Simpson  = simpson_index
) %>%
  left_join(
    data.frame(sample_data(ps_norm)),
    by = "SampleID"
  )

write.csv(richness_diversity_manual, "richness_diversity.csv",
          quote = FALSE, row.names = FALSE)

richness_diversity_manual <- richness_diversity_manual %>%
  mutate(across(c(Sample, rep, Lesion, Host, Location,
                  Group, Group_lesion, Sample_lesion), factor))

## 7a. Shannon diversity (linear model) ####
model_shannon <- lm(Shannon ~ Group_lesion, data = richness_diversity_manual)

res_shannon <- DHARMa::simulateResiduals(model_shannon) |>
  resid(quantileFunction = qnorm, outlierValues = c(-5, 5))
ggResidpanel::resid_auxpanel(res_shannon, fitted(model_shannon))

emmeans::joint_tests(model_shannon)
emm_shannon <- emmeans::emmeans(model_shannon, ~ Group_lesion)
cld_results_sha  <- multcomp::cld(emm_shannon, Letters = letters, adjust = "none")
cld_results_sha$.group <- trimws(cld_results_sha$.group)
p_shannon <- cld_results_sha |>
  as.data.frame() |>
  ggplot(aes(x = Group_lesion, fill = Group_lesion)) +
  geom_col(aes(y = emmean)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  geom_text(aes(y = emmean + SE + 0.1, label = .group), size = 5) +
  theme_Publication() +
  labs(x = "", y = "Shannon diversity", fill = "Lesion-host group") +
  scale_fill_manual(values = host_colours) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_shannon
## 7b. Observed richness (negative binomial GLM) ####
poisson_model   <- glm(Observed ~ Group_lesion, family = poisson,
                       data = richness_diversity_manual)
model_observed1 <- glmmTMB(Observed ~ Group_lesion,
                           family = nbinom1(link = "log"),
                           data = richness_diversity_manual)
model_observed2 <- glmmTMB(Observed ~ Group_lesion,
                           family = nbinom2(link = "log"),
                           data = richness_diversity_manual)
AIC(poisson_model, model_observed1, model_observed2)

model_observed <- model_observed1   # nbinom1 selected on lowest AIC

res_observed <- DHARMa::simulateResiduals(model_observed) |>
  resid(quantileFunction = qnorm, outlierValues = c(-5, 5))
ggResidpanel::resid_auxpanel(res_observed, fitted(model_observed))

emmeans::joint_tests(model_observed)
emm_observed <- emmeans::emmeans(model_observed, ~ Group_lesion, type = "response")
cld_results  <- multcomp::cld(emm_observed, Letters = letters, adjust = "none")
cld_results$.group <- trimws(cld_results$.group)

p_observed <- as.data.frame(cld_results) |>
  ggplot(aes(x = Group_lesion, fill = Group_lesion)) +
  geom_col(aes(y = response)) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE), width = 0.2) +
  geom_text(aes(y = response + SE + 5, label = .group), size = 5) +
  theme_Publication() +
  labs(x = "", y = "Observed richness", fill = "Lesion-host group") +
  scale_fill_manual(values = host_colours) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


# 8. Beta diversity ####

## 8a. Bray-Curtis (on log-transformed counts) ####
logt        <- transform_sample_counts(ps_norm, function(x) log(1 + x))
out.pcoa.bc <- ordinate(logt, method = "MDS", distance = "bray")

p_pcoa_bc <- plot_ordination(logt, out.pcoa.bc,
                             color = "Group_lesion") +
  geom_point(
    size = 2,
    position = position_jitter(width = 0.05, height = 0.05),
    alpha = 0.8
  ) +
  theme_Publication() +
  scale_color_manual(values = host_colours)

ps_bray      <- phyloseq::distance(logt, method = "bray")
meta_beta <- data.frame(sample_data(logt))
set.seed(123)
(ado_t_bray  <- vegan::adonis2(ps_bray ~ Group_lesion, data = meta_beta))
set.seed(123)
(ado_l_bray  <- vegan::adonis2(ps_bray ~ Host, data = meta_beta))
set.seed(123)
(ado_tl_bray <- vegan::adonis2(ps_bray ~ Host / Group_lesion, data = meta_beta))
set.seed(123)
beta_disp_bc <- betadisper(ps_bray, sample_data(logt)$Group_lesion)
set.seed(123)
anova(beta_disp_bc)
set.seed(123)
permutest(beta_disp_bc, permutations = 999)

## 8b. Unweighted UniFrac ####

out.pcoa.uu <- ordinate(
  ps_norm,
  method = "MDS",
  distance = "uunifrac"
)

p_pcoa_uu <- plot_ordination(
  ps_norm,
  out.pcoa.uu,
  color = "Group_lesion"
) +
  geom_point(
    size = 2,
    position = position_jitter(width = 0.05, height = 0.05),
    alpha = 0.8
  ) +
  theme_Publication() +
  scale_color_manual(values = host_colours)

ps_uunifrac <- phyloseq::distance(ps_norm, method = "uunifrac")

meta_uu <- data.frame(sample_data(ps_norm))

set.seed(123)
(ado_t_uu <- vegan::adonis2(
  ps_uunifrac ~ Group_lesion,
  data = meta_uu))
set.seed(123)
(ado_l_uu <- vegan::adonis2(
  ps_uunifrac ~ Host,
  data = meta_uu))
set.seed(123)
(ado_tl_uu <- vegan::adonis2(
  ps_uunifrac ~ Host / Group_lesion,
  data = meta_uu))
set.seed(123)
(beta_disp_uu <- vegan::betadisper(
  ps_uunifrac,
  meta_uu$Group_lesion
))
set.seed(123)
anova(beta_disp_uu)
set.seed(123)
permutest(beta_disp_uu, permutations = 999)

## 8c. Assemble diversity figure ####

p_pcoa_bc <- p_pcoa_bc +
  labs(title = "Bray-Curtis")

p_pcoa_uu <- p_pcoa_uu +
  labs(title = "Unweighted UniFrac")

diversity_plots <- ggarrange(
  p_observed, p_shannon, p_pcoa_bc, p_pcoa_uu,
  nrow = 2, ncol = 2,
  common.legend = TRUE, legend = "right",
  align = "hv", labels = c("A", "", "B", "")
)
diversity_plots
ggsave("Fig2_diversity_plots.png", diversity_plots, width = 8, height = 5, bg = "white")


# 9. Community dispersion analysis ####
# Quantifies within-group variability (distance to group centroid) using
# the Bray-Curtis betadisper objects already computed in Section 8.
# Complements the PERMANOVA results above by testing whether differences
# in community composition (Section 8) reflect true compositional shifts
# versus differences in within-group heterogeneity.

## 9a. Dispersion by lesion-host group ####
anova_beta_disp_group     <- anova(beta_disp_bc)
permutest_beta_disp_group <- permutest(beta_disp_bc, permutations = 999)

print(anova_beta_disp_group)
print(permutest_beta_disp_group)

disp_group_df <- data.frame(
  SampleID            = names(beta_disp_bc$distances),
  DistanceToCentroid  = beta_disp_bc$distances
) %>%
  left_join(
    data.frame(sample_data(ps_norm)) %>%
      rownames_to_column("SampleID_check") %>%
      select(-SampleID_check),  # or just confirm SampleID already exists and drop the rownames version
    by = "SampleID"
  )

write_csv(disp_group_df, "bray_distance_to_centroid_group_lesion.csv")

disp_group_summary <- disp_group_df %>%
  group_by(Group_lesion, Host) %>%
  summarise(
    mean_distance = mean(DistanceToCentroid),
    se_distance   = sd(DistanceToCentroid) / sqrt(n()),
    n             = n(),
    .groups = "drop"
  )

write_csv(disp_group_summary, "bray_distance_to_centroid_summary.csv")
print(disp_group_summary)

p_disp_group <- ggplot(
  disp_group_df,
  aes(x = Group_lesion, y = DistanceToCentroid, fill = Group_lesion)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = host_colours) +
  theme_bw() +
  labs(x = "", y = "Bray-Curtis distance to centroid",
       title = "Community dispersion among lesion groups") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave("community_dispersion_group_lesion.png", p_disp_group,
       width = 6, height = 4, bg = "white")

## 9b. Dispersion by host ####
beta_disp_host <- betadisper(ps_bray, sample_data(logt)$Host)

anova_beta_disp_host     <- anova(beta_disp_host)
permutest_beta_disp_host <- permutest(beta_disp_host, permutations = 999)

print(anova_beta_disp_host)
print(permutest_beta_disp_host)

disp_host_df <- data.frame(
  SampleID           = names(beta_disp_host$distances),
  DistanceToCentroid = beta_disp_host$distances
) %>%
  left_join(
    as.data.frame(sample_data(logt)),
    by = "SampleID"
  )
write_csv(disp_host_df, "bray_distance_to_centroid_host.csv")

p_disp_host <- ggplot(
  disp_host_df,
  aes(x = Host, y = DistanceToCentroid, fill = Host)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  theme_bw() +
  labs(x = "", y = "Bray-Curtis distance to centroid",
       title = "Host-level community dispersion") +
  theme(legend.position = "none")

ggsave("community_dispersion_host.png", p_disp_host,
       width = 5, height = 4, bg = "white")

#10. Sampling type####
ps_gs      <- subset_samples(ps_norm, Host == "Goosegrass")
logt_gs    <- transform_sample_counts(ps_gs, function(x) log(1 + x))
ps_bray_gs <- phyloseq::distance(logt_gs, method = "bray")
meta_gs    <- as(sample_data(ps_gs), "data.frame")

out.pcoa.gs <- ordinate(logt_gs, method = "MDS", distance = "bray")
set.seed(123)
ado_lesion_gs <- vegan::adonis2(
  ps_bray_gs ~ Location + Lesion, data = meta_gs, by = "margin"
)

# Test whether the effect of sampling strategy differs by location
ado_lesion_interaction_gs <- vegan::adonis2(
  ps_bray_gs ~ Location * Lesion,
  data = meta_gs,
  permutations = 999,
  by = "margin"
)

print(ado_lesion_interaction_gs)

cat("\n=== PERMANOVA: Lesion type within goosegrass ===\n")
print(ado_lesion_gs)

beta_disp_lesion <- vegan::betadisper(ps_bray_gs, meta_gs$Lesion)
cat("\n=== PERMDISP: Lesion type within goosegrass ===\n")
print(anova(beta_disp_lesion))
print(permutest(beta_disp_lesion, permutations = 999))

p_pcoa_lesion <- plot_ordination(logt_gs, out.pcoa.gs,
                                 color = "Lesion") +
  geom_point(aes(shape = Location), size = 3) +
  stat_ellipse(aes(group = Lesion, color = Lesion),
               type = "t", linetype = 2, alpha = 0.4) +
  scale_color_manual(values = c("PL" = "#DC267F", "SL" = "#FE6100")) +
  labs(color = "Lesion type", shape = "Location") +
  theme_Publication()

ggsave("pcoa_goosegrass_lesion_type.png",
       p_pcoa_lesion, width = 6, height = 5, bg = "white")
# 11. Core microbiome analysis ####
# Defines "core" as ASVs present in at least a specified proportion of samples
# within each host and above a minimum relative abundance threshold.

core_prevalence_threshold <- 0.50   # present in >=50% of samples within each host
core_abundance_threshold  <- 0.001  # >=0.1% relative abundance in a sample

# Relative abundance table from raw filtered counts
ps_rel_core <- transform_sample_counts(ps1, function(x) x / sum(x))

otu_rel_core <- as(otu_table(ps_rel_core), "matrix")

if (!taxa_are_rows(ps_rel_core)) {
  otu_rel_core <- t(otu_rel_core)
}

storage.mode(otu_rel_core) <- "numeric"

meta_core <- data.frame(sample_data(ps_rel_core))
meta_core$SampleID <- rownames(meta_core)

stopifnot(all(colnames(otu_rel_core) %in% meta_core$SampleID))

host_levels <- unique(meta_core$Host)

# Presence is defined as relative abundance >= threshold
otu_core_presence <- otu_rel_core >= core_abundance_threshold

core_prevalence_by_host <- lapply(host_levels, function(h) {
  
  samples_h <- meta_core$SampleID[meta_core$Host == h]
  
  data.frame(
    ASV = rownames(otu_core_presence),
    Host = h,
    Prevalence = rowMeans(otu_core_presence[, samples_h, drop = FALSE]),
    MeanRelativeAbundance = rowMeans(otu_rel_core[, samples_h, drop = FALSE]),
    MaxRelativeAbundance = apply(
      otu_rel_core[, samples_h, drop = FALSE],
      1,
      max
    )
  )
}) %>%
  bind_rows()

write_csv(
  core_prevalence_by_host,
  "core_prevalence_by_host.csv"
)

# Host-specific core ASVs
core_by_host <- core_prevalence_by_host %>%
  filter(Prevalence >= core_prevalence_threshold) %>%
  group_by(Host) %>%
  summarise(
    Core_ASVs = list(ASV),
    n_core = n(),
    .groups = "drop"
  )

write_csv(
  core_by_host %>% select(Host, n_core),
  "core_counts_by_host.csv"
)

# Shared core = core in BOTH hosts
core_lists <- core_by_host$Core_ASVs
names(core_lists) <- core_by_host$Host

shared_core_asvs <- Reduce(intersect, core_lists)

# Host-specific cores
goosegrass_core_only <- setdiff(
  core_lists[["Goosegrass"]],
  core_lists[["Corn"]]
)

corn_core_only <- setdiff(
  core_lists[["Corn"]],
  core_lists[["Goosegrass"]]
)

core_taxonomy <- data.frame(
  ASV = unique(c(shared_core_asvs, goosegrass_core_only, corn_core_only))
) %>%
  mutate(
    CoreCategory = case_when(
      ASV %in% shared_core_asvs ~ "Shared core",
      ASV %in% goosegrass_core_only ~ "Goosegrass core only",
      ASV %in% corn_core_only ~ "Corn core only",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(
    as.data.frame(tax_table(ps1)) %>%
      rownames_to_column("ASV"),
    by = "ASV"
  ) %>%
  left_join(
    core_prevalence_by_host %>%
      select(ASV, Host, Prevalence, MeanRelativeAbundance) %>%
      pivot_wider(
        names_from = Host,
        values_from = c(Prevalence, MeanRelativeAbundance)
      ),
    by = "ASV"
  ) %>%
  arrange(CoreCategory, ASV)

write_csv(core_taxonomy, "core_microbiome_taxonomy.csv")

cat("\nCore microbiome definition:\n")
cat("Prevalence threshold:", core_prevalence_threshold, "\n")
cat("Abundance threshold:", core_abundance_threshold, "\n")

cat("\nShared core ASVs:", length(shared_core_asvs), "\n")
cat("Corn core-only ASVs:", length(corn_core_only), "\n")
cat("Goosegrass core-only ASVs:", length(goosegrass_core_only), "\n")

print(core_taxonomy)

# Per-sample relative abundance of each core category
core_abundance_by_sample <- data.frame(
  SampleID = colnames(otu_rel_core),
  SharedCore_RA = colSums(
    otu_rel_core[shared_core_asvs, , drop = FALSE],
    na.rm = TRUE
  ),
  CornCoreOnly_RA = colSums(
    otu_rel_core[corn_core_only, , drop = FALSE],
    na.rm = TRUE
  ),
  GoosegrassCoreOnly_RA = colSums(
    otu_rel_core[goosegrass_core_only, , drop = FALSE],
    na.rm = TRUE
  )
) %>%
  left_join(meta_core, by = "SampleID")

write_csv(
  core_abundance_by_sample,
  "core_relative_abundance_by_sample.csv"
)

core_abundance_long <- core_abundance_by_sample %>%
  pivot_longer(
    cols = c(SharedCore_RA, CornCoreOnly_RA, GoosegrassCoreOnly_RA),
    names_to = "CoreCategory",
    values_to = "RelativeAbundance"
  ) %>%
  mutate(
    CoreCategory = recode(
      CoreCategory,
      SharedCore_RA = "Shared core",
      CornCoreOnly_RA = "Corn core only",
      GoosegrassCoreOnly_RA = "Goosegrass core only"
    )
  )

core_abundance_summary <- core_abundance_long %>%
  group_by(Host, Group_lesion, CoreCategory) %>%
  summarise(
    MeanRelativeAbundance = mean(RelativeAbundance),
    SERelativeAbundance = sd(RelativeAbundance) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

write_csv(
  core_abundance_summary,
  "core_relative_abundance_summary.csv"
)

p_core_abundance <- ggplot(
  core_abundance_long,
  aes(x = Group_lesion, y = RelativeAbundance, fill = Group_lesion)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  facet_wrap(~ CoreCategory, scales = "free_y") +
  scale_fill_manual(values = host_colours) +
  theme_bw() +
  labs(
    x = "",
    y = "Relative abundance",
    title = "Core microbiome relative abundance"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  "core_microbiome_relative_abundance.png",
  p_core_abundance,
  width = 8,
  height = 4,
  bg = "white"
)

## Figure 4A ####
otu_pa_fig4 <- as(otu_table(ps1), "matrix")

if (!taxa_are_rows(ps1)) {
  otu_pa_fig4 <- t(otu_pa_fig4)
}

otu_pa_fig4 <- otu_pa_fig4 > 0

meta_fig4 <- data.frame(sample_data(ps1))
meta_fig4$SampleID <- rownames(meta_fig4)

# Identify ASVs present within each lesion-host group

asv_sets <- lapply(unique(meta_fig4$Group_lesion), function(g) {
  
  samples_g <- meta_fig4$SampleID[
    meta_fig4$Group_lesion == g
  ]
  
  rownames(otu_pa_fig4)[
    rowSums(
      otu_pa_fig4[, samples_g, drop = FALSE]
    ) > 0
  ]
})

names(asv_sets) <- unique(meta_fig4$Group_lesion)

asv_sets

# Create presence/absence table for UpSet plot

all_asvs <- unique(unlist(asv_sets))

upset_df <- data.frame(
  ASV = all_asvs
)

for (grp in names(asv_sets)) {
  
  upset_df[[grp]] <- all_asvs %in% asv_sets[[grp]]
  
}

head(upset_df)

# Desired group order

group_order <- c(
  "Corn PL",
  "Corn SL",
  "Goosegrass PL",
  "Goosegrass SL"
)

# Convert to list-column format required by ggupset

upset_long <- upset_df %>%
  pivot_longer(
    cols = all_of(group_order),
    names_to = "Group_lesion",
    values_to = "Present"
  ) %>%
  filter(Present) %>%
  group_by(ASV) %>%
  summarise(
    Groups = list(Group_lesion),
    .groups = "drop"
  )

# Figure 4A

p_upset_simple <- ggplot(
  upset_long,
  aes(x = Groups)
) +
  
  # Intersection-size bars
  geom_bar(
    fill = "grey40",
    width = 0.85
  ) +
  
  # Add intersection count above each bar
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.35,
    size = 3.5
  ) +
  
  # Order intersections from largest to smallest
  ggupset::scale_x_upset(
    order_by = "freq"
  ) +
  
  # Give labels enough room above tallest bar
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.10)
    )
  ) +
  
  labs(
    x = "Lesion-host group intersection",
    y = "Number of ASVs"
  ) +
  
  theme_bw() +
  theme_Publication()

p_upset_simple

## Figure 4B: Core + non-core stacked mean relative abundance ####

# Colors matched to biological categories
core_cols <- c(
  "Corn core only"       = "#648FFF",
  "Shared core"          = "#785EF0",
  "Goosegrass core only" = "#FE6100",
  "Non-core ASVs"        = "grey75"
)

# Desired stacking order: bottom to top
core_stack_order <- c(
  "Corn core only",
  "Shared core",
  "Goosegrass core only",
  "Non-core ASVs"
)

# Desired legend order
core_legend_order <- c(
  "Corn core only",
  "Shared core",
  "Goosegrass core only",
  "Non-core ASVs"
)

core_abundance_summary_plot <- core_abundance_long %>%
  mutate(
    Group_lesion = factor(
      Group_lesion,
      levels = c("Corn PL", "Corn SL", "Goosegrass PL", "Goosegrass SL")
    ),
    CoreCategory = factor(
      CoreCategory,
      levels = c("Corn core only", "Shared core", "Goosegrass core only")
    )
  ) %>%
  group_by(SampleID, Group_lesion, CoreCategory) %>%
  summarise(
    RelativeAbundance = sum(RelativeAbundance),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = CoreCategory,
    values_from = RelativeAbundance,
    values_fill = 0
  ) %>%
  mutate(
    `Non-core ASVs` =
      1 - (`Corn core only` + `Shared core` + `Goosegrass core only`),
    `Non-core ASVs` = pmax(`Non-core ASVs`, 0)
  ) %>%
  pivot_longer(
    cols = all_of(core_stack_order),
    names_to = "CoreCategory",
    values_to = "RelativeAbundance"
  ) %>%
  group_by(Group_lesion, CoreCategory) %>%
  summarise(
    MeanRelativeAbundance = mean(RelativeAbundance),
    .groups = "drop"
  ) %>%
  mutate(
    CoreCategory = factor(CoreCategory, levels = core_stack_order)
  )

p_core_stacked <- ggplot(
  core_abundance_summary_plot,
  aes(
    x = Group_lesion,
    y = MeanRelativeAbundance,
    fill = CoreCategory
  )
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = core_cols,
    breaks = core_legend_order,
    name = NULL
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = NULL,
    y = "Mean relative abundance"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(t = 0, unit = "cm"),
    legend.box.spacing = unit(0.1, "cm"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, color = "black"),
    plot.margin = margin(5, 5, 0, 5)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             override.aes = list(color = "black", linewidth = 0.25)))

p_core_stacked

## Figure 4C: Host-level indicator ASVs ####

fig4c_df <- ind_host_sig %>%
  mutate(
    Host_association = case_when(
      index == 1 ~ "Corn",
      index == 2 ~ "Goosegrass",
      TRUE ~ NA_character_
    ),
    
    Taxon_label = case_when(
      !is.na(Species) & Species != "" ~
        paste0(
          "<i>", Genus, " ", Species, "</i>",
          " (", ASV, ")"
        ),
      TRUE ~
        paste0(
          "<i>", Genus, "</i>",
          " (", ASV, ")"
        )
    ),
    
    Signed_IndVal = ifelse(
      Host_association == "Corn",
      -stat,
      stat
    ),
    
    Stat_label = paste0(
      sprintf("%.2f", stat),
      "\n(p = ",
      sprintf("%.3f", p.value),
      ")"
    ),
    
    # Put label immediately outside its corresponding point
    Label_x = ifelse(
      Host_association == "Corn",
      Signed_IndVal - 0.035,
      Signed_IndVal + 0.035
    ),
    
    # Justify whole two-line label away from point
    Label_hjust = ifelse(
      Host_association == "Corn",
      1,
      0
    )
  ) %>%
  arrange(Signed_IndVal) %>%
  mutate(
    Taxon_label = factor(
      Taxon_label,
      levels = Taxon_label
    )
  )


fig4c <- ggplot(
  fig4c_df,
  aes(
    x = Signed_IndVal,
    y = Taxon_label
  )
) +
  
  # Zero reference
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey55"
  ) +
  
  # Lollipop stems
  geom_segment(
    aes(
      x = 0,
      xend = Signed_IndVal,
      yend = Taxon_label,
      color = Host_association
    ),
    linewidth = 1.1
  ) +
  
  # Points
  geom_point(
    aes(color = Host_association),
    size = 4
  ) +
  
  # Statistics: use explicit x coordinate
  geom_text(
    data = fig4c_df,
    aes(
      x = Label_x,
      y = Taxon_label,
      label = Stat_label,
      hjust = Label_hjust
    ),
    inherit.aes = FALSE,
    size = 3.3,
    lineheight = 0.95,
    color = "black"
  ) +
  
  # Host headings
  annotate(
    "text",
    x = -0.5,
    y = Inf,
    label = "Corn",
    color = host_colours_host["Corn"],
    fontface = "bold",
    vjust = 0.3,
    size = 4.2
  ) +
  
  annotate(
    "text",
    x = 0.5,
    y = Inf,
    label = "Goosegrass",
    color = host_colours_host["Goosegrass"],
    fontface = "bold",
    vjust = 0.3,
    size = 4.2
  ) +
  
  scale_color_manual(
    values = host_colours_host
  ) +
  
  # Give statistics room outside ±1
  scale_x_continuous(
    limits = c(-1.32, 1.32),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = c("1.0", "0.5", "0", "0.5", "1.0"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  labs(
    x = "Indicator value (IndVal)",
    y = NULL
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    # No redundant legend
    legend.position = "none",
    
    axis.text.y = ggtext::element_markdown(
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 11
    ),
    
    # Give host labels space above panel
    plot.margin = margin(
      t = 15,
      r = 10,
      b = 5,
      l = 5
    )
  )

fig4c

fig4 <- (
  p_upset_simple + p_core_stacked +
    plot_layout(widths = c(1.35, 1))
) /
  wrap_elements(fig4c) +
  plot_layout(
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

fig4

ggsave(
  "Figure4_shared_ASVs_core_microbiome.png",
  fig4,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure4_shared_ASVs_core_microbiome.svg",
  fig4,
  width = 11,
  height = 9
)

# 12. Phylogenetic tree of top 50 ASVs with bootstrap support####

select <- dplyr::select

# 1. Select top 50 ASVs by total abundance

otu_ps1 <- as(otu_table(ps1), "matrix")

if (!taxa_are_rows(ps1)) {
  otu_ps1 <- t(otu_ps1)
}

top50 <- rowSums(otu_ps1) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(50)

ps50 <- prune_taxa(top50, ps1)


# 2. Build top-50 reference tree with bootstrap support

tree50 <- ape::keep.tip(
  fitGTR$tree,
  top50
)

bs50 <- lapply(
  bs,
  function(tr) {
    ape::keep.tip(tr, top50)
  }
)

class(bs50) <- "multiPhylo"

cat("Number of bootstrap trees:", length(bs50), "\n")

# Recalculate bootstrap support for the exact top-50 tree
tree50_boot <- plotBS(
  tree50,
  bs50,
  type = "none"
)

cat("Bootstrap labels:\n")
print(tree50_boot$node.label)

# Save tree with bootstrap support
ape::write.tree(
  tree50_boot,
  file = "top50_ASV_tree_bootstrap_newick.tre"
)

# 3. Build ggtree object with bootstrap labels

p_base <- ggtree(
  tree50_boot,
  ladderize = TRUE
) +
  geom_text2(
    aes(
      subset =
        !isTip &
        !is.na(suppressWarnings(as.numeric(label))) &
        suppressWarnings(as.numeric(label)) >= 0.70,
      label = round(
        suppressWarnings(as.numeric(label)) * 100,
        0
      )
    ),
    hjust = 1.2,
    vjust = -0.35,
    size = 2.5
  )

tree_edges <- p_base$data

tree_tips <- tree_edges %>%
  filter(isTip) %>%
  transmute(
    ASV = label,
    x_tip = x,
    y_tip = y
  )

# 4. Taxonomy labels

tax50_clean <- as.data.frame(tax_table(ps50)) %>%
  rownames_to_column("ASV") %>%
  mutate(
    Species_clean = case_when(
      is.na(Species) |
        Species == "" |
        Species == "Uncharacterized" ~ "sp.",
      TRUE ~ Species
    ),
    
    Genus_clean = case_when(
      is.na(Genus) |
        Genus == "" |
        Genus == "Uncharacterized" ~ "Uncharacterized",
      TRUE ~ Genus
    ),
    
    label_plot = case_when(
      Species_clean == "sp." ~
        paste0(
          "'",
          ASV,
          " '~italic('",
          Genus_clean,
          "')~' sp.'"
        ),
      
      TRUE ~
        paste0(
          "'",
          ASV,
          " '~italic('",
          Genus_clean,
          "')~italic('",
          Species_clean,
          "')"
        )
    )
  ) %>%
  left_join(
    tree_tips,
    by = "ASV"
  )

# 5. Presence/absence data for top 50 ASVs

otu50 <- as(otu_table(ps50), "matrix")

if (!taxa_are_rows(ps50)) {
  otu50 <- t(otu50)
}

meta50 <- data.frame(sample_data(ps50))
meta50$SampleID <- rownames(meta50)

group_order <- c(
  "Corn PL",
  "Corn SL",
  "Goosegrass PL",
  "Goosegrass SL"
)

stopifnot(
  all(group_order %in% meta50$Group_lesion)
)

presence_long <- lapply(
  seq_along(group_order),
  function(i) {
    
    g <- group_order[i]
    
    samples_g <- meta50$SampleID[
      meta50$Group_lesion == g
    ]
    
    data.frame(
      ASV = rownames(otu50),
      Group_lesion = g,
      group_i = i,
      Present = rowSums(
        otu50[, samples_g, drop = FALSE]
      ) > 0
    )
  }
) %>%
  bind_rows() %>%
  left_join(
    tree_tips,
    by = "ASV"
  )

# 6. Classify each ASV by host association

presence_class <- presence_long %>%
  mutate(
    Host = case_when(
      grepl("^Corn", Group_lesion) ~ "Corn",
      grepl("^Goosegrass", Group_lesion) ~ "Goosegrass"
    )
  ) %>%
  group_by(ASV) %>%
  summarise(
    present_corn =
      any(Present & Host == "Corn"),
    
    present_goosegrass =
      any(Present & Host == "Goosegrass"),
    
    .groups = "drop"
  ) %>%
  mutate(
    Presence_category = case_when(
      present_corn &
        !present_goosegrass ~ "Corn only",
      
      present_corn &
        present_goosegrass ~ "Shared",
      
      !present_corn &
        present_goosegrass ~ "Goosegrass only",
      
      TRUE ~ "Not detected"
    )
  )

presence_long <- presence_long %>%
  left_join(
    presence_class %>%
      select(
        ASV,
        Presence_category
      ),
    by = "ASV"
  ) %>%
  mutate(
    Heatmap_fill = if_else(
      Present,
      Presence_category,
      "Absent"
    )
  )


# 7. Heatmap and label coordinates

label_x_offset <- 0.01

heatmap_cell_width  <- 0.058
heatmap_cell_height <- 0.78

heatmap_x_start <-
  max(tree_tips$x_tip, na.rm = TRUE) + 0.5

presence_long <- presence_long %>%
  mutate(
    heatmap_x =
      heatmap_x_start +
      (group_i - 1) * heatmap_cell_width
  )

# Labels directly above individual heatmap columns
# Labels directly above individual heatmap columns
heatmap_labels <- data.frame(
  label = c("PL", "SL", "PL", "SL"),
  heatmap_x =
    heatmap_x_start +
    (seq_along(group_order) - 1) * heatmap_cell_width,
  y_lab =
    max(tree_tips$y_tip, na.rm = TRUE) + 0.55
)

# Host labels centered over each pair of columns
host_labels <- data.frame(
  label = c("C", "G"),
  heatmap_x = c(
    mean(
      heatmap_x_start +
        c(0, 1) * heatmap_cell_width
    ),
    mean(
      heatmap_x_start +
        c(2, 3) * heatmap_cell_width
    )
  ),
  y_lab =
    max(tree_tips$y_tip, na.rm = TRUE) + 1.8
)

# Horizontal segments spanning each pair of heatmap cells
segment_inset <- heatmap_cell_width * 0.15

host_segments <- data.frame(
  x = c(
    heatmap_x_start - heatmap_cell_width / 2 + segment_inset,
    heatmap_x_start + 2 * heatmap_cell_width - heatmap_cell_width / 2 + segment_inset
  ),
  xend = c(
    heatmap_x_start + heatmap_cell_width + heatmap_cell_width / 2 - segment_inset,
    heatmap_x_start + 3 * heatmap_cell_width + heatmap_cell_width / 2 - segment_inset
  ),
  y = rep(
    max(tree_tips$y_tip, na.rm = TRUE) + 1.6,
    2
  ),
  yend = rep(
    max(tree_tips$y_tip, na.rm = TRUE) + 1.6,
    2
  )
)

# Sample abbreviation legend
sample_legend_df <- data.frame(
  Sample_key = factor(
    c("C = Corn",
      "G = Goosegrass",
      "PL = Pooled lesion",
      "SL = Single lesion"),
    levels = c(
      "C = Corn",
      "G = Goosegrass",
      "PL = Pooled lesion",
      "SL = Single lesion"
    )
  ),
  x = heatmap_x_start,
  y = max(tree_tips$y_tip, na.rm = TRUE) + 20
)
# 8. Family colors

families <- unique(tax50_clean$Family)

family_cols <- scales::hue_pal()(
  length(families)
)

names(family_cols) <- families

family_legend_df <- tax50_clean %>%
  distinct(Family) %>%
  mutate(
    x = heatmap_x_start,
    y = max(tree_tips$y_tip, na.rm = TRUE) + 20
  )

# 9. Host-association palette

presence_cols <- c(
  "Absent"          = "grey96",
  "Corn only"       = "#648FFF",
  "Shared"          = "#785EF0",
  "Goosegrass only" = "#FE6100"
)

# 10. Final tree + colored heatmap + bootstrap values

family_breaks <- names(family_cols)

p_tree_heat <- p_base +
  
  # Dummy points for Family legend
  geom_point(
    data = family_legend_df,
    aes(
      x = x,
      y = y,
      color = Family
    ),
    shape = 16,
    size = 3,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  
  # Taxonomy labels
  geom_text(
    data = tax50_clean,
    aes(
      x = x_tip + label_x_offset,
      y = y_tip,
      label = label_plot,
      color = Family
    ),
    parse = TRUE,
    hjust = 0,
    size = 3.9,
    show.legend = FALSE
  ) +
  
  # Heatmap
  geom_tile(
    data = presence_long,
    aes(
      x = heatmap_x,
      y = y_tip,
      fill = Heatmap_fill
    ),
    width = heatmap_cell_width,
    height = heatmap_cell_height,
    color = "grey85",
    linewidth = 0.18
  ) +
  
  # Horizontal host grouping segments
  geom_segment(
    data = host_segments,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.7,
    color = "black"
  ) +
  
  # PL / SL labels
  geom_text(
    data = heatmap_labels,
    aes(
      x = heatmap_x,
      y = y_lab,
      label = label
    ),
    hjust = 0.5,
    vjust = 0,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  
  # C / G labels
  geom_text(
    data = host_labels,
    aes(
      x = heatmap_x,
      y = y_lab,
      label = label
    ),
    hjust = 0.5,
    vjust = 0,
    size = 4,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    values = family_cols,
    breaks = family_breaks,
    name = "Family"
  ) +
  
  scale_fill_manual(
    values = presence_cols,
    breaks = c(
      "Absent",
      "Corn only",
      "Shared",
      "Goosegrass only"
    ),
    name = "Host association"
  ) +
  
  coord_cartesian(
    xlim = c(
      0,
      max(presence_long$heatmap_x, na.rm = TRUE)
    ),
    ylim = c(
      min(tree_tips$y_tip, na.rm = TRUE) - 1,
      max(tree_tips$y_tip, na.rm = TRUE) + 2.5
    ),
    clip = "off"
  ) +
  
  theme_tree2() +
  
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.justification = "left",
    
    text = element_text(size = 11),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11, face = "bold"),
    
    plot.margin = margin(12, 25, 5, 5)
  ) +
  
  guides(
    fill = guide_legend(
      title = "Host association",
      order = 1,
      override.aes = list(
        color = "grey80"
      )
    ),
    
    color = guide_legend(
      title = "Family",
      order = 2,
      override.aes = list(
        shape = 16,
        size = 3,
        alpha = 1,
        colour = unname(family_cols[family_breaks])
      )
    )
  )


library(cowplot)

#
# Sample key
#

sample_key <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 1,
    label = paste0(
      "Sample\n",
      "C = Corn\n",
      "G = Goosegrass\n",
      "PL = Pooled lesion\n",
      "SL = Single lesion"
    ),
    hjust = 0,
    vjust = 1,
    size = 3.7,
    lineheight = 1.15
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )


#
# Extract Host association + Family legends
#

legend <- cowplot::get_legend(
  p_tree_heat +
    theme(
      legend.position = "right",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.spacing.y = grid::unit(2, "pt")
    )
)


#
# Remove legends from main tree
#

tree_only <- p_tree_heat +
  theme(
    legend.position = "none"
  )


#
# Build right-hand panel manually
#
# This avoids the large blank space produced by plot_grid()
#

right_panel <- ggdraw() +
  
  # Sample block at the top
  draw_plot(
    sample_key,
    x = 0.00,
    y = 0.80,
    width = 1.00,
    height = 0.18
  ) +
  
  # Host association + Family immediately below
  draw_plot(
    legend,
    x = 0.00,
    y = 0.20,
    width = 1.00,
    height = 0.62
  )


#
# Final figure
#

p_tree_final <- plot_grid(
  tree_only,
  right_panel,
  nrow = 1,
  rel_widths = c(
    0.80,
    0.20
  ),
  align = "h"
)

p_tree_final

ggsave(
  "Fig3_top50_tree_bootstrap_host_heatmap.png",
  p_tree_final,
  width = 320,
  height = 220,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  "Fig3_top50_tree_bootstrap_host_heatmap.svg",
  p_tree_final,
  width = 320,
  height = 220,
  units = "mm",
  dpi = 600,
  bg = "white"
)

### Abudance ####
otu_mat  <- as.matrix(otu_table(ps50))
tax_df   <- as.data.frame(tax_table(ps50)) %>% rownames_to_column("ASV")
asv_sums <- colSums(otu_mat)
tax_abundance <- tax_df %>%
  mutate(TotalAbundance = asv_sums[ASV]) %>%
  group_by(Family) %>%
  summarise(TotalAbundance = sum(TotalAbundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  mutate(PercentReads = TotalAbundance / sum(TotalAbundance) * 100)


# 13. Differential abundance analysis (edgeR quasi-likelihood F-test) ####

otu_matrix_da <- as(otu_table(ps1), "matrix")
if (!taxa_are_rows(ps1)) otu_matrix_da <- t(otu_matrix_da)

group <- factor(sample_data(ps1)$Host, levels = c("Corn", "Goosegrass"))

y_da <- DGEList(counts = otu_matrix_da, group = group)
y_da <- calcNormFactors(y_da)

design <- model.matrix(~ group)
y_da <- estimateDisp(y_da, design)
fit_da <- glmQLFit(y_da, design)
qlf <- glmQLFTest(fit_da)

degs <- topTags(qlf, n = Inf)$table %>%
  rownames_to_column("ASV") %>%
  left_join(
    as.data.frame(tax_table_original) %>% rownames_to_column("ASV"),
    by = "ASV"
  )

sig_asvs <- degs %>% filter(FDR < 0.05)

write_csv(degs,     "degs.csv")
write_csv(sig_asvs, "sig_asvs.csv")

sig_asvs_clean <- sig_asvs %>%
  mutate(
    `Log2 Fold Change`       = round(logFC, 2),
    `Adjusted p-value (FDR)` = signif(FDR, 3),
    Taxonomy                 = paste(Genus, Species)
  ) %>%
  select(ASV, `Log2 Fold Change`, `Adjusted p-value (FDR)`, Taxonomy)

write_csv(sig_asvs_clean, "sig_asvs_clean.csv")


# 14. Indicator species analysis ####
# Complements the differential abundance analysis (Section 13) using the
# IndVal approach (Dufrene and Legendre 1997), which identifies taxa
# significantly associated with a particular group based on both
# specificity and fidelity, rather than fold-change in abundance alone.

## 14a. Indicator species by host ####
# Relative abundance
ps_rel_ind <- transform_sample_counts(ps1, function(x) x / sum(x))

# OTU/ASV table: samples as rows
otu_ind <- as(otu_table(ps_rel_ind), "matrix")

if (taxa_are_rows(ps_rel_ind)) {
  otu_ind <- t(otu_ind)
}

otu_ind <- as.data.frame(otu_ind)

# Host grouping variable, matched to OTU table sample order
meta_ind <- data.frame(sample_data(ps_rel_ind))

group_ind <- meta_ind[rownames(otu_ind), "Host"]
group_ind <- factor(group_ind)

# Check
table(group_ind)
identical(rownames(otu_ind), sample_names(ps_rel_ind))

# Indicator species analysis
set.seed(123)

ind_host <- multipatt(
  otu_ind,
  group_ind,
  func = "IndVal.g",
  control = how(nperm = 999)
)

ind_host_df <- ind_host$sign %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  left_join(
    as.data.frame(tax_table(ps_norm)) %>% rownames_to_column("ASV"),
    by = "ASV"
  ) %>%
  arrange(p.value)

write_csv(ind_host_df, "indicator_species_host.csv")

ind_host_sig <- ind_host_df %>% filter(p.value <= 0.05)
write_csv(ind_host_sig, "indicator_species_host_significant.csv")
print(ind_host_sig)

## 14b. Indicator species by lesion-host group ####
group_ind_lesion <- sample_data(ps_norm)$Group_lesion

set.seed(123)
ind_group_lesion <- multipatt(
  otu_ind, group_ind_lesion,
  func    = "IndVal.g",
  control = how(nperm = 999)
)
summary(ind_group_lesion)

ind_group_lesion_df <- ind_group_lesion$sign %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  left_join(
    as.data.frame(tax_table(ps_norm)) %>% rownames_to_column("ASV"),
    by = "ASV"
  ) %>%
  arrange(p.value)

write_csv(ind_group_lesion_df, "indicator_species_group_lesion.csv")

ind_group_lesion_sig <- ind_group_lesion_df %>% filter(p.value <= 0.05)
write_csv(ind_group_lesion_sig, "indicator_species_group_lesion_significant.csv")
print(ind_group_lesion_sig)


# 15. Phyllachora maydis identity analysis ####

# Identify ASVs sharing >= 97% sequence identity with ASV3 (P. maydis)
ASV3_seq   <- as.character(ASVs.nochim["ASV3"])
pid_values <- sapply(ASVs.nochim, function(s) {
  pwalign::pid(pwalign::pairwiseAlignment(pattern = ASV3_seq, subject = s))
})
high_id_ASVs <- names(which(pid_values >= 97))
cat("ASVs >= 97% identity with ASV3:", paste(high_id_ASVs, collapse = ", "), "\n")

seqs_Pmaydis <- ASVs.nochim[high_id_ASVs]
writeXStringSet(seqs_Pmaydis, "pmaydis.fasta", format = "fasta")

selected_ASVs <- t(tmp.seqtab[, colnames(tmp.seqtab) %in% high_id_ASVs])
host_groups   <- split(colnames(selected_ASVs), samdf$Host)
total_reads_pm <- data.frame(
  ASV        = rownames(selected_ASVs),
  Corn       = rowSums(selected_ASVs[, host_groups$Corn,       drop = FALSE]),
  Goosegrass = rowSums(selected_ASVs[, host_groups$Goosegrass, drop = FALSE])
)
write.csv(total_reads_pm, "pmaydis_asvs_reads.csv",
          row.names = FALSE, quote = FALSE)

alignment_pm <- msa(seqs_Pmaydis)
aligned_pm   <- msaConvert(alignment_pm, type = "ape::DNAbin")
pi_obs       <- nuc.div(aligned_pm)

set.seed(123)
n_permutations <- 1000
bootstrap_pi   <- replicate(n_permutations, {
  idx <- sample(seq_len(nrow(aligned_pm)), replace = TRUE)
  nuc.div(aligned_pm[idx, ])
})
p_value_pi <- mean(bootstrap_pi >= pi_obs)
cat("Observed nucleotide diversity (pi):", round(pi_obs, 4), "\n")
cat("Bootstrap p-value:", round(p_value_pi, 3), "\n")


# P. luzhouensis ASVs based on >= 97% identity to ASV2

# Identify ASVs sharing >= 97% sequence identity with ASV2 (P. luzhouensis)
ASV2_seq <- as.character(ASVs.nochim["ASV2"])

pid_values_pl <- sapply(ASVs.nochim, function(s) {
  pwalign::pid(
    pwalign::pairwiseAlignment(
      pattern = ASV2_seq,
      subject = s
    )
  )
})

high_id_ASVs_pl <- names(which(pid_values_pl >= 97))

cat(
  "ASVs >= 97% identity with ASV2:",
  paste(high_id_ASVs_pl, collapse = ", "),
  "\n"
)



# Export sequences


seqs_Pluzhouensis <- ASVs.nochim[high_id_ASVs_pl]

writeXStringSet(
  seqs_Pluzhouensis,
  "pluzhouensis.fasta",
  format = "fasta"
)


# Read counts by host

selected_ASVs_pl <- t(
  tmp.seqtab[, colnames(tmp.seqtab) %in% high_id_ASVs_pl]
)

host_groups_pl <- split(
  colnames(selected_ASVs_pl),
  samdf$Host
)

total_reads_pl <- data.frame(
  ASV = rownames(selected_ASVs_pl),
  
  Corn = rowSums(
    selected_ASVs_pl[, host_groups_pl$Corn, drop = FALSE]
  ),
  
  Goosegrass = rowSums(
    selected_ASVs_pl[, host_groups_pl$Goosegrass, drop = FALSE]
  )
)

write.csv(
  total_reads_pl,
  "pluzhouensis_asvs_reads.csv",
  row.names = FALSE,
  quote = FALSE
)


# Multiple sequence alignment


alignment_pl <- msa(seqs_Pluzhouensis)

aligned_pl <- msaConvert(
  alignment_pl,
  type = "ape::DNAbin"
)


# Nucleotide diversity


pi_obs_pl <- nuc.div(aligned_pl)


# Bootstrap analysis


set.seed(123)

n_permutations <- 1000

bootstrap_pi_pl <- replicate(
  n_permutations,
  {
    idx <- sample(
      seq_len(nrow(aligned_pl)),
      replace = TRUE
    )
    
    nuc.div(aligned_pl[idx, ])
  }
)

p_value_pi_pl <- mean(
  bootstrap_pi_pl >= pi_obs_pl
)


# Results

cat(
  "Observed nucleotide diversity (pi) for P. luzhouensis:",
  round(pi_obs_pl, 4),
  "\n"
)

cat(
  "Bootstrap p-value:",
  round(p_value_pi_pl, 3),
  "\n"
)

cat(
  "Number of ASVs >= 97% identity with ASV2:",
  length(high_id_ASVs_pl),
  "\n"
)

print(total_reads_pl)





#
# RESULTS OUTPUT FOR MANUSCRIPT SECTION 3.2
# Fungal Community Associated with Tar Spot Lesions
#

select <- dplyr::select

cat("\n",
    "============================================================\n",
    "3.2 FUNGAL COMMUNITY ASSOCIATED WITH TAR SPOT LESIONS\n",
    "============================================================\n",
    sep = "")


#
# 1. Dataset summary
#

cat("\n=== 1. DATASET SUMMARY ===\n")

cat("Samples in filtered dataset:", nsamples(ps1), "\n")
cat("ASVs retained after filtering:", ntaxa(ps1), "\n")

cat("\nSamples by host and lesion type:\n")
print(
  data.frame(sample_data(ps1)) %>%
    count(Host, Group_lesion)
)


#
# 2. Alpha diversity
#

cat("\n=== 2. ALPHA DIVERSITY ===\n")

cat("\nObserved richness - joint test:\n")
print(emmeans::joint_tests(model_observed))

cat("\nObserved richness - estimated marginal means:\n")
print(as.data.frame(emm_observed))

cat("\nShannon diversity - joint test:\n")
print(emmeans::joint_tests(model_shannon))

cat("\nShannon diversity - estimated marginal means:\n")
print(as.data.frame(emm_shannon))


#
# 3. Beta diversity: Bray-Curtis
#

cat("\n=== 3. BETA DIVERSITY: BRAY-CURTIS ===\n")

cat("\nPERMANOVA - lesion-host group:\n")
print(ado_t_bray)

cat("\nPERMANOVA - host:\n")
print(ado_l_bray)

cat("\nPERMANOVA - host / lesion-host group nested model:\n")
print(ado_tl_bray)


#
# 4. Beta diversity: unweighted UniFrac
#

cat("\n=== 4. BETA DIVERSITY: UNWEIGHTED UNIFRAC ===\n")

cat("\nPERMANOVA - lesion-host group:\n")
print(ado_t_uu)

cat("\nPERMANOVA - host:\n")
print(ado_l_uu)

cat("\nPERMANOVA - host / lesion-host group nested model:\n")
print(ado_tl_uu)


#
# 5. Community dispersion
#

cat("\n=== 5. COMMUNITY DISPERSION: BRAY-CURTIS ===\n")

cat("\nPERMDISP - lesion-host group, ANOVA:\n")
print(anova_beta_disp_group)

cat("\nPERMDISP - lesion-host group, permutation test:\n")
print(permutest_beta_disp_group)

cat("\nMean distance to centroid by lesion-host group:\n")
print(disp_group_summary)

cat("\nPERMDISP - host, ANOVA:\n")
print(anova_beta_disp_host)

cat("\nPERMDISP - host, permutation test:\n")
print(permutest_beta_disp_host)


#
# 6. Taxonomic composition of the top 50 ASVs
#

cat("\n=== 6. TAXONOMIC COMPOSITION OF TOP 50 ASVs ===\n")

# Recalculate robustly with ASVs as rows
otu50_results <- as(otu_table(ps50), "matrix")

if (!taxa_are_rows(ps50)) {
  otu50_results <- t(otu50_results)
}

tax50_results <- as.data.frame(tax_table(ps50)) %>%
  rownames_to_column("ASV") %>%
  mutate(
    TotalReads = rowSums(otu50_results)[ASV]
  )

family_abundance_results <- tax50_results %>%
  group_by(Family) %>%
  summarise(
    TotalReads = sum(TotalReads),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalReads)) %>%
  mutate(
    PercentReads = 100 * TotalReads / sum(TotalReads)
  )

print(family_abundance_results)


#
# 7. Phyllachora ASVs assigned directly by taxonomy
#

cat("\n=== 7. TAXONOMICALLY ASSIGNED PHYLLACHORA ASVs ===\n")

phyllachora_results <- as.data.frame(tax_table(ps1)) %>%
  rownames_to_column("ASV") %>%
  filter(Genus == "Phyllachora") %>%
  select(ASV, Family, Genus, Species)

print(phyllachora_results)

cat(
  "\nNumber of ASVs directly assigned to Phyllachora:",
  nrow(phyllachora_results),
  "\n"
)


#
# 8. Host distribution of ASV2 and ASV3
#

cat("\n=== 8. DOMINANT PHYLLACHORA ASVs BY HOST ===\n")

otu_host_results <- as(otu_table(ps1), "matrix")

if (!taxa_are_rows(ps1)) {
  otu_host_results <- t(otu_host_results)
}

meta_host_results <- data.frame(sample_data(ps1))
meta_host_results$SampleID <- rownames(meta_host_results)

host_read_summary <- lapply(
  c("ASV2", "ASV3"),
  function(a) {
    
    data.frame(
      ASV = a,
      Corn = sum(
        otu_host_results[
          a,
          meta_host_results$SampleID[
            meta_host_results$Host == "Corn"
          ],
          drop = FALSE
        ]
      ),
      Goosegrass = sum(
        otu_host_results[
          a,
          meta_host_results$SampleID[
            meta_host_results$Host == "Goosegrass"
          ],
          drop = FALSE
        ]
      )
    )
  }
) %>%
  bind_rows() %>%
  left_join(
    as.data.frame(tax_table(ps1)) %>%
      rownames_to_column("ASV") %>%
      select(ASV, Genus, Species),
    by = "ASV"
  )

print(host_read_summary)


#
# 9. >=97% identity groups around ASV3 and ASV2
#

cat("\n=== 9. PHYLLACHORA SEQUENCE GROUPS (>=97% ITS1 IDENTITY) ===\n")

cat("\nP. maydis-like ASVs (>=97% identity to ASV3):\n")
cat(paste(high_id_ASVs, collapse = ", "), "\n")

cat("Number of P. maydis-like ASVs:", length(high_id_ASVs), "\n")

cat("\nRead counts by host:\n")
print(total_reads_pm)

cat(
  "Total P. maydis-like reads in corn:",
  sum(total_reads_pm$Corn),
  "\n"
)

cat(
  "Total P. maydis-like reads in goosegrass:",
  sum(total_reads_pm$Goosegrass),
  "\n"
)

cat(
  "Nucleotide diversity (pi), P. maydis-like group:",
  round(pi_obs, 4),
  "\n"
)


cat("\nP. luzhouensis-like ASVs (>=97% identity to ASV2):\n")
cat(paste(high_id_ASVs_pl, collapse = ", "), "\n")

cat(
  "Number of P. luzhouensis-like ASVs:",
  length(high_id_ASVs_pl),
  "\n"
)

cat("\nRead counts by host:\n")
print(total_reads_pl)

cat(
  "Total P. luzhouensis-like reads in corn:",
  sum(total_reads_pl$Corn),
  "\n"
)

cat(
  "Total P. luzhouensis-like reads in goosegrass:",
  sum(total_reads_pl$Goosegrass),
  "\n"
)

cat(
  "Nucleotide diversity (pi), P. luzhouensis-like group:",
  round(pi_obs_pl, 4),
  "\n"
)


#
# 10. Indicator species by host
#

cat("\n=== 10. INDICATOR SPECIES BY HOST ===\n")

ind_host_results <- ind_host_sig %>%
  select(
    ASV,
    Genus,
    Species,
    Family,
    stat,
    p.value,
    index
  )

print(ind_host_results)

cat(
  "\nNumber of significant host indicator ASVs:",
  nrow(ind_host_results),
  "\n"
)


#
# 11. Indicator species by lesion-host group
#

cat("\n=== 11. INDICATOR SPECIES BY LESION-HOST GROUP ===\n")

ind_group_results <- ind_group_lesion_sig %>%
  select(
    ASV,
    Genus,
    Species,
    Family,
    stat,
    p.value,
    index
  )

print(ind_group_results)

cat(
  "\nNumber of significant lesion-host indicator ASVs:",
  nrow(ind_group_results),
  "\n"
)


#
# 12. ASV overlap among lesion-host groups
#

cat("\n=== 12. ASV OVERLAP AMONG LESION-HOST GROUPS ===\n")

cat("\nNumber of ASVs detected in each lesion-host group:\n")
print(
  tibble(
    Group_lesion = names(asv_sets),
    n_ASVs = lengths(asv_sets)
  )
)

shared_all_four <- Reduce(intersect, asv_sets)

cat(
  "\nASVs detected across all four lesion-host groups:",
  length(shared_all_four),
  "\n"
)

cat("Identity of ASVs detected across all four groups:\n")
print(shared_all_four)

shared_all_four_tax <- as.data.frame(tax_table(ps1)) %>%
  rownames_to_column("ASV") %>%
  filter(ASV %in% shared_all_four) %>%
  select(ASV, Family, Genus, Species)

print(shared_all_four_tax)


# Host-level presence/absence overlap
corn_samples_overlap <- meta_fig4$SampleID[
  meta_fig4$Host == "Corn"
]

goose_samples_overlap <- meta_fig4$SampleID[
  meta_fig4$Host == "Goosegrass"
]

corn_detected_asvs <- rownames(otu_pa_fig4)[
  rowSums(
    otu_pa_fig4[, corn_samples_overlap, drop = FALSE]
  ) > 0
]

goose_detected_asvs <- rownames(otu_pa_fig4)[
  rowSums(
    otu_pa_fig4[, goose_samples_overlap, drop = FALSE]
  ) > 0
]

shared_host_asvs <- intersect(
  corn_detected_asvs,
  goose_detected_asvs
)

cat(
  "\nASVs detected in both corn and goosegrass:",
  length(shared_host_asvs),
  "\n"
)


#
# 13. Core microbiome
#

cat("\n=== 13. HOST-LEVEL CORE MICROBIOME ===\n")

cat(
  "Core definition: relative abundance >=",
  core_abundance_threshold * 100,
  "% in at least",
  core_prevalence_threshold * 100,
  "% of samples within a host\n"
)

cat(
  "\nShared core ASVs:",
  length(shared_core_asvs),
  "\n"
)

cat(
  "Corn core-only ASVs:",
  length(corn_core_only),
  "\n"
)

cat(
  "Goosegrass core-only ASVs:",
  length(goosegrass_core_only),
  "\n"
)

cat("\nCore ASV taxonomy:\n")

print(
  core_taxonomy %>%
    select(
      ASV,
      CoreCategory,
      Family,
      Genus,
      Species,
      starts_with("Prevalence"),
      starts_with("MeanRelativeAbundance")
    )
)

cat("\nMean relative abundance of core categories:\n")
print(core_abundance_summary)


#
# 14. Top-50 Phyllachora host-association pattern
#

cat("\n=== 14. PHYLLACHORA IN THE TOP-50 PHYLOGENY / HEATMAP ===\n")

top50_phyllachora <- tax50_clean %>%
  filter(Genus_clean == "Phyllachora") %>%
  select(ASV, Family, Genus_clean, Species_clean) %>%
  left_join(
    presence_class %>%
      select(
        ASV,
        present_corn,
        present_goosegrass,
        Presence_category
      ),
    by = "ASV"
  )

print(top50_phyllachora)

cat(
  "\nNumber of Phyllachora ASVs represented among top 50:",
  nrow(top50_phyllachora),
  "\n"
)


cat("\n",
    "============================================================\n",
    "END OF SECTION 3.2 OUTPUT\n",
    "============================================================\n",
    sep = "")

