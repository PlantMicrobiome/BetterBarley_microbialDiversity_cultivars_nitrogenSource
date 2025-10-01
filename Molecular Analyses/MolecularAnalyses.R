############################################################
# Microbiome analysis pipeline (ITS / protists / 16S)
# - Loads OTU & taxonomy tables
# - Cleans taxonomy strings and filters data
# - Builds phyloseq objects
# - Rarefies 16S data and aggregates to Genus/Phylum
# - Creates relative-abundance bar plot (top 10 genera + Others)
# - Trains RF to identify important genera and plots heatmap
# - Computes alpha diversity (Observed, Chao1, Shannon, Simpson)
# - Fits simple linear models vs. inoculum
# - Computes beta diversity (Bray–Curtis) and PCoA + adonis2
#
# Notes:
# * Expects QIIME-style taxonomy strings like "d__Bacteria;p__Firmicutes;..."
# * Expects metadata to include columns: substrate, plant, ino (and plant_sub for PCoA)
# * Set your working directory so all relative paths resolve.
############################################################

# ---------------------------- #
# 0) Libraries
# ---------------------------- #

# Core DE/phylogenetics/plots/utilities
library("DESeq2");      packageVersion("DESeq2")
library("ape");         packageVersion("ape")
library("phyloseq");    packageVersion("phyloseq")
library("ggplot2");     packageVersion("ggplot2")
library("decontam");    packageVersion("decontam")
library("tidyverse");   packageVersion("tidyverse")
library("dplyr");       packageVersion("dplyr")
library("vegan");       packageVersion("vegan")

# Microbiome helpers / palettes / heatmaps
library("MicrobiotaProcess")
library("microbiome")
library(Polychrome)
library(tidyr)
library(pheatmap)

# Named palette for genera (RDS must contain named vector: names == genera)
palette_named <- readRDS("genus_palette.rds")
colorblind_Palette = c("#DDAA33", "#4DDFD3", "#004488", "#FF6600", "#AA4499")
# ---------------------------- #
# 1) Load feature tables (OTU/ASV)
# ---------------------------- #
# Each table: rows = taxa (OTUs/ASVs), columns = samples; row 1 is header.
otu_its  <- read.csv("its_OTU.csv",         row.names = 1, sep = "\t", header = TRUE)
otu_pro  <- read.csv("protist_OTU.csv",     row.names = 1, sep = "\t", header = TRUE)
otu_16s  <- read.csv("16s_OTU.csv",         row.names = 1, sep = "\t", header = TRUE)

# ---------------------------- #
# 2) Load taxonomy tables
# ---------------------------- #
# Expect a 'Taxon' column with semicolon-delimited ranks (d__/p__/c__/...).
taxonomy_its <- read.table("its_taxonomy.tsv",     sep = "\t", header = TRUE, row.names = 1)
taxonomy_pro <- read.table("protist_taxonomy.tsv", sep = "\t", header = TRUE, row.names = 1)
taxonomy_16s <- read.table("16s_taxonomy.tsv",     sep = "\t", header = TRUE, row.names = 1)

# ---------------------------- #
# 3) Clean taxonomy strings -> (Kingdom..Species) columns
# ---------------------------- #
# Strips “d__”, “p__”, etc. and splits into 7 ranks. Missing levels retained as NA/""
# Tip: If some rows have fewer levels, 'separate()' will fill with NAs.
# ---- ITS ----
tax_its <- taxonomy_its %>%
  select(Taxon) %>%
  separate(Taxon, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), ";", fill = "right")

tax_its.clean <- data.frame(
  row.names = row.names(tax_its),
  Kingdom = str_replace(tax_its[,1], "d__", ""),
  Phylum  = str_replace(tax_its[,2], "p__", ""),
  Class   = str_replace(tax_its[,3], "c__", ""),
  Order   = str_replace(tax_its[,4], "o__", ""),
  Family  = str_replace(tax_its[,5], "f__", ""),
  Genus   = str_replace(tax_its[,6], "g__", ""),
  Species = str_replace(tax_its[,7], "s__", ""),
  stringsAsFactors = FALSE
)

# ---- Protists ----
tax_pro <- taxonomy_pro %>%
  select(Taxon) %>%
  separate(Taxon, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), ";", fill = "right")

tax_pro.clean <- data.frame(
  row.names = row.names(tax_pro),
  Kingdom = str_replace(tax_pro[,1], "d__", ""),
  Phylum  = str_replace(tax_pro[,2], "p__", ""),
  Class   = str_replace(tax_pro[,3], "c__", ""),
  Order   = str_replace(tax_pro[,4], "o__", ""),
  Family  = str_replace(tax_pro[,5], "f__", ""),
  Genus   = str_replace(tax_pro[,6], "g__", ""),
  Species = str_replace(tax_pro[,7], "s__", ""),
  stringsAsFactors = FALSE
)

# ---- 16S ----
tax_16s <- taxonomy_16s %>%
  select(Taxon) %>%
  separate(Taxon, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), ";", fill = "right")

tax_16s.clean <- data.frame(
  row.names = row.names(tax_16s),
  Kingdom = str_replace(tax_16s[,1], "d__", ""),
  Phylum  = str_replace(tax_16s[,2], "p__", ""),
  Class   = str_replace(tax_16s[,3], "c__", ""),
  Order   = str_replace(tax_16s[,4], "o__", ""),
  Family  = str_replace(tax_16s[,5], "f__", ""),
  Genus   = str_replace(tax_16s[,6], "g__", ""),
  Species = str_replace(tax_16s[,7], "s__", ""),
  stringsAsFactors = FALSE
)

# ---------------------------- #
# 4) Load metadata (QIIME format)
# ---------------------------- #
# Must contain rownames matching sample IDs in OTU tables.
mapfile_its <- import_qiime_sample_data("its_metadata.tsv")
mapfile_pro <- import_qiime_sample_data("protist_metadata.tsv")
mapfile_16s <- import_qiime_sample_data("16s_metadata.tsv")

# ---------------------------- #
# 5) Build phyloseq objects
# ---------------------------- #
OTU_ITS <- otu_table(as.matrix(otu_its), taxa_are_rows = TRUE)
OTU_PRO <- otu_table(as.matrix(otu_pro), taxa_are_rows = TRUE)
OTU_16S <- otu_table(as.matrix(otu_16s), taxa_are_rows = TRUE)

TAX_ITS <- phyloseq::tax_table(as.matrix(tax_its.clean))
TAX_PRO <- phyloseq::tax_table(as.matrix(tax_pro.clean))
TAX_16S <- phyloseq::tax_table(as.matrix(tax_16s.clean))

META_ITS <- sample_data(mapfile_its)
META_PRO <- sample_data(mapfile_pro)
META_16S <- sample_data(mapfile_16s)

physeq_its  <- phyloseq(OTU_ITS,  TAX_ITS,  META_ITS)
physeq_pro  <- phyloseq(OTU_PRO,  TAX_PRO,  META_PRO)
physeq_16s  <- phyloseq(OTU_16S,  TAX_16S,  META_16S)

# ---------------------------- #
# 6) Taxonomic filtering / decontamination heuristics
# ---------------------------- #
# ITS: remove non-Fungi, uncultured/Incertaesedis/empty/unwanted plant algae, etc.
clean_its <- subset_taxa(physeq_its, Phylum != "NA")
clean_its <- subset_taxa(clean_its, Genus  != "NA")
clean_its <- subset_taxa(clean_its, Genus  != "uncultured")
clean_its <- subset_taxa(clean_its, Phylum != "Fungi_phy_Incertae_sedis")
clean_its <- subset_taxa(clean_its, Genus  != "Fungi_gen_Incertae_sedis")
clean_its <- subset_taxa(clean_its, Genus  != "")
clean_its <- subset_taxa(clean_its, Phylum != "Streptophyta") # plant
clean_its <- subset_taxa(clean_its, Phylum != "Chlorophyta")  # algae
clean_its <- subset_taxa(clean_its, Genus  != "Sporidiobolaceae_gen_Incertae_sedis")
clean_its <- subset_taxa(clean_its, Genus  != "Lentitheciaceae_gen_Incertae_sedis")
clean_its <- subset_taxa(clean_its, Genus  != "GS11_gen_Incertae_sedis")
clean_its <- subset_taxa(clean_its, Genus  != "Malasseziaceae_gen_Incertae_sedis")
clean_its <- subset_taxa(clean_its, Genus  != "Sebacinales_gen_Incertae_sedis")
clean_its <- subset_taxa(clean_its, Genus  != "Candida")
clean_its <- subset_taxa(clean_its, (Kingdom == "Fungi") | is.na(Kingdom))

# Prevalence filter: keep taxa present (>1 read) in >0.5% of samples
final_clean_its <- filter_taxa(clean_its, function(x) sum(x > 1) > (0.005 * length(x)), prune = TRUE)
final_clean_its

# Protists: keep Cercozoa; remove NA/uncultured
clean_pro <- subset_taxa(physeq_pro, Phylum != "NA")
clean_pro <- subset_taxa(clean_pro, Genus  != "NA")
clean_pro <- subset_taxa(clean_pro, Phylum == "Cercozoa")
clean_pro <- subset_taxa(clean_pro, Genus  != "uncultured")
final_clean_pro <- filter_taxa(clean_pro, function(x) sum(x > 1) > (0.005 * length(x)), prune = TRUE)

# 16S: drop chloroplast/mitochondria; remove NA/uncultured; prevalence filter
clean_16s <- subset_taxa(physeq_16s, (Order  != "Chloroplast")  | is.na(Order))
clean_16s <- subset_taxa(clean_16s, (Family != "Mitochondria") | is.na(Family))
clean_16s <- subset_taxa(clean_16s, Genus  != "NA")
clean_16s <- subset_taxa(clean_16s, Phylum != "NA")
clean_16s <- subset_taxa(clean_16s, Genus  != "uncultured")
final_clean_16s <- filter_taxa(clean_16s, function(x) sum(x > 1) > (0.005 * length(x)), prune = TRUE)

# ---------------------------- #
# 7) 16S rarefaction and agglomeration
# ---------------------------- #
# Rarefy to even depth for downstream ordination/plotting comparability.
# NOTE: 'rngseed' expects an integer; TRUE will be coerced to 1. Prefer a fixed seed (e.g., 123).
bac_rare <- rarefy_even_depth(final_clean_16s, rngseed = 123)

# Agglomerate counts to Phylum and Genus (keeps NA/"" if NArm=FALSE)
bac_rare_phylum <- tax_glom(bac_rare, taxrank = "Phylum", NArm = FALSE, bad_empty = c(NA, "", " ", "\t"))
bac_rare_genus  <- tax_glom(bac_rare, taxrank = "Genus",  NArm = FALSE, bad_empty = c(NA, "", " ", "\t"))

# ---------------------------- #
# 8) Figure 3 — Relative abundance barplot (top 10 genera per treatment)
# ---------------------------- #
# Compute % relative abundance per sample, then average by (substrate, plant, ino).
bac_rel <- transform_sample_counts(bac_rare, function(x) x / sum(x) * 100)

# Aggregate to genus and melt for tidy operations
bac_rel_genus      <- tax_glom(bac_rel, taxrank = "Genus")
bac_rel_genus_melt <- psmelt(bac_rel_genus)

# Mean abundance per group (Treatment = plant_ino) to stabilize across replicates
bac_rel_genus_melt_sum <- bac_rel_genus_melt %>%
  group_by(substrate, plant, ino, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = 'drop') %>%
  mutate(Treatment = paste(plant, ino, sep = "_"))

# Identify top 10 genera within each (substrate, plant, ino) combination
top10_genera <- bac_rel_genus_melt_sum %>%
  group_by(substrate, plant, ino) %>%
  slice_max(order_by = Abundance, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Top10 = TRUE)

# Label non-top10 as "Others" and collapse them
melt_top10 <- bac_rel_genus_melt_sum %>%
  left_join(select(top10_genera, substrate, plant, ino, Genus, Top10),
            by = c("substrate", "plant", "ino", "Genus")) %>%
  mutate(Genus = ifelse(is.na(Top10), "Others", as.character(Genus))) %>%
  group_by(substrate, plant, ino, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(Treatment = paste(plant, ino, sep = "_"))

# Genus factor levels: alphabetical with "Others" placed last
all_genera   <- unique(melt_top10$Genus) %>% as.character()
genus_levels <- c(sort(setdiff(all_genera, "Others")), "Others")
melt_top10$Genus <- factor(melt_top10$Genus, levels = genus_levels)

# Stacked bar plot faceted by substrate (rows) and plant (columns)
barplot_top10 <- ggplot(melt_top10, aes(x = ino, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_grid(substrate ~ plant, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = palette_named) +
  labs(x = "ino", y = "Mean Relative Abundance (%)") +
  theme_minimal() +
  theme(
    legend.position = "left",
    axis.text.y     = element_text(size = 16, color = "black"),
    axis.title.x    = element_blank(),
    axis.text.x     = element_blank(),
    axis.title.y    = element_text(size = 16, color = "black", face = "bold"),
    strip.text      = element_text(size = 16, color = "black"),
    panel.grid      = element_blank(),
    panel.background= element_rect(fill = "white", colour = "white"),
    axis.line.y     = element_line(color = "black", linewidth = 0.8),
    axis.ticks.y    = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length.y = unit(0.3, "cm")
  ) +
  labs(y = "Relative Abundance (%)")

barplot_top10
ggsave('barplot_top10.png', barplot_top10, width = 12, height = 10, dpi = 768)

# ---------------------------- #
# 9) Figure 4 — Random Forest: important genera + heatmap
# ---------------------------- #
library(randomForest)

# Create wide matrix: one row per sample, one column per Genus, values = relative abundance
abundance_wide <- bac_rel_genus_melt %>%
  ungroup() %>%
  select(sample = Sample, Genus, Abundance) %>%
  pivot_wider(names_from = Genus, values_from = Abundance, values_fill = 0)

# Ensure valid column names for modeling
colnames(abundance_wide) <- make.names(colnames(abundance_wide))

# Treatment label per sample (plant_ino then joined with substrate later)
sample_treatment <- bac_rel_genus_melt %>%
  ungroup() %>%
  mutate(Treatment = paste(substrate, plant, ino, sep = "_")) %>%
  select(sample = Sample, Treatment) %>%
  distinct()

rf_data <- left_join(sample_treatment, abundance_wide, by = "sample")

# Train RF to classify Treatment using genera abundances
set.seed(123)
rf_model <- randomForest(as.factor(Treatment) ~ .,
                         data = rf_data %>% select(-sample),
                         ntree = 50000, importance = TRUE)

# Extract feature importance and select top 15 by MeanDecreaseAccuracy
imp_genera     <- importance(rf_model)
imp_genera_df  <- data.frame(Genus = rownames(imp_genera),
                             MeanDecreaseGini = imp_genera[, "MeanDecreaseGini"],
                             MeanDecreaseAccuracy = imp_genera[, "MeanDecreaseAccuracy"])
top_15_important <- imp_genera_df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  slice_head(n = 15)

# Subset to top genera and average abundance across identical Treatments
top_genera_names <- top_15_important$Genus

abundance_top <- rf_data %>%
  select(sample, Treatment, all_of(top_genera_names))

abundance_avg <- abundance_top %>%
  group_by(Treatment) %>%
  summarise_at(vars(all_of(top_genera_names)), mean) %>%
  ungroup()

# Matrix with Treatment as rows; z-score per genus (row-scale)
abundance_avg_mat <- abundance_avg %>%
  column_to_rownames(var = "Treatment") %>%
  as.matrix()

# Column annotations derived from Treatment string (substrate, plant, ino)
annotation_df <- data.frame(Treatment = rownames(abundance_avg_mat)) %>%
  mutate(
    substrate = sapply(strsplit(as.character(Treatment), "_"), `[`, 1),
    plant     = sapply(strsplit(as.character(Treatment), "_"), `[`, 2),
    ino       = sapply(strsplit(as.character(Treatment), "_"), `[`, 3)
  ) %>%
  select(-Treatment) %>%
  mutate(across(everything(), ~factor(.)))

rownames(annotation_df) <- rownames(abundance_avg_mat)

# Z-score by genus (rows)
abundance_scaled <- t(scale(t(abundance_avg_mat)))

# Annotation color maps (adjust to your palette needs)
substrate_colors <- c(Chitin = "#b892ff", Mineral = "#4cc9f0")
plant_colors     <- c(Babushka = "#80ed99", RGT = "#ffafcc")
ino_colors       <- c("10-1" = "#ff686b", "10-3" = "#ff9b85", "10-5" = "#ffd97d",
                      "10-7" = "#aaf683", "no. ino" = "#60d394")

annotation_colors <- list(substrate = substrate_colors,
                          plant     = plant_colors,
                          ino       = ino_colors)

# Order annotation columns for legend display
annotation_df <- annotation_df %>% select(ino, plant, substrate)

# Heatmap: columns = Treatments, rows = genera (after transpose)
heatmap_transposed <- pheatmap(
  t(abundance_scaled),
  annotation_col   = annotation_df,
  cluster_rows     = FALSE,
  cluster_cols     = FALSE,
  show_rownames    = TRUE,
  show_colnames    = FALSE,
  annotation_colors = annotation_colors
)

ggsave('heatmap_zscore_unclustered.png', heatmap_transposed, width = 6, height = 3, dpi = 600)

# ---------------------------- #
# 10) Figure 5 — Alpha diversity vs. inoculum
# ---------------------------- #
# Compute richness metrics on genus-level rarefied table
alpha_bac_rare <- estimate_richness(bac_rare_genus,
                                    measures = c("Observed", "Shannon", "Chao1", "Simpson"))

# Design (sample_data) extracted from phyloseq object
design <- bac_rare_genus@sam_data

# Extract design columns by position (consider switching to name-based selection to avoid brittle indexing)
design_substrate <- as.data.frame(design[, 6]); rownames(design_substrate) <- rownames(design); colnames(design_substrate) <- "Substrate"
design_cultivar  <- as.data.frame(design[, 4]); rownames(design_cultivar)  <- rownames(design); colnames(design_cultivar)  <- "Cultivar"
design_inoculum  <- as.data.frame(design[, 5]); rownames(design_inoculum)  <- rownames(design); colnames(design_inoculum)  <- "Inoculum"

designFinal <- cbind(design_cultivar, design_inoculum, design_substrate)

# ---------------- Observed richness ----------------
alpha_bac_rare_Observed <- as.data.frame(alpha_bac_rare[ , "Observed"])
colnames(alpha_bac_rare_Observed) <- "Observed"

alpha_bac_rare_Observed_final <- cbind(designFinal, alpha_bac_rare) %>% as.data.frame()

# Order factors for consistent plotting
alpha_bac_rare_Observed_final$Substrate <- ordered(alpha_bac_rare_Observed_final$Substrate, levels = c("Chitin","Mineral"))
alpha_bac_rare_Observed_final$Inoculum <- ordered(alpha_bac_rare_Observed_final$Inoculum, levels = c("1","3","5","7","9"))
alpha_bac_rare_Observed_final$Cultivar <- ordered(alpha_bac_rare_Observed_final$Cultivar, levels = c("Babushka","RGT"))

# rename maps for facets/legends if needed
cultivar_names  <- c(Babushka = "Babushka", RGT = "RGT")
substrate_names <- c(Chitin = "Chitin", Mineral = "Mineral")
inoculum_names  <- c(`10-1`="1", `10-3`="3", `10-5`="5", `10-7`="7", `no. ino`="9")

# ---------------- Chao1 ----------------
alpha_bac_rare_Chao1 <- as.data.frame(alpha_bac_rare[ , "Chao1"])
colnames(alpha_bac_rare_Chao1) <- "Chao1"
alpha_bac_rare_Chao1_final <- cbind(designFinal, alpha_bac_rare_Chao1) %>% as.data.frame()

alpha_bac_rare_Chao1_final$Substrate <- ordered(alpha_bac_rare_Chao1_final$Substrate, levels = c("Chitin","Mineral"))
alpha_bac_rare_Chao1_final$Inoculum <- ordered(alpha_bac_rare_Chao1_final$Inoculum, levels = c("1","3","5","7","9"))
alpha_bac_rare_Chao1_final$Cultivar <- ordered(alpha_bac_rare_Chao1_final$Cultivar, levels = c("Babushka","RGT"))

# Axis labels for alpha-diversity plots
x.axis.label <- "Microbiome"
y.axis.label <- "Chao1"

# Convert ordered factors to numeric positions (1,3,5,7,9) for regression lines
alpha_bac_rare_Chao1_final$Inoculum <- as.numeric(levels(alpha_bac_rare_Chao1_final$Inoculum))[alpha_bac_rare_Chao1_final$Inoculum]

# Plot Chao1 across inoculum dilutions, faceted by substrate
chao1_plot <- ggplot(alpha_bac_rare_Chao1_final, aes(Inoculum, Chao1, color = Cultivar, shape = Cultivar)) +
  # Two-layer smooth lines per cultivar (bold with SE, then thin without SE for edge effect)
  geom_smooth(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "Babushka"), color = "#ff791f", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "RGT"),      color = "#000000", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "Babushka"), color = "#000000", method = "lm", se = FALSE, linewidth = 1) +
  geom_smooth(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "RGT"),      color = "#FFFFFF", method = "lm", se = FALSE, linewidth = 1) +
  # Points with outline
  geom_point(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "Babushka"), shape = 16, size = 4, color = "#ff791f") +
  geom_point(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "RGT"),      shape = 17, size = 4, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "Babushka"), shape = 16, size = 3, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "RGT"),      shape = 17, size = 3, color = "#FFFFFF") +
  labs(x = x.axis.label, y = y.axis.label) +
  facet_grid(~ Substrate) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.justification = "right",
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", face = "bold"),
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    legend.text   = element_text(size = 12, color = "black"),
    legend.title  = element_blank(),
    strip.text.x  = element_text(size = 12, color = "black"),
    strip.text.y  = element_text(size = 12, color = "black"),
    panel.spacing = unit(1.5, "lines")
  ) +
  scale_x_continuous(
    breaks = c(1,3,5,7,9),
    labels = c(expression(10^-1), expression(10^-3), expression(10^-5), expression(10^-7), expression("No"))
  ) +
  scale_shape_manual(values = c("Babushka" = 16, "RGT" = 17))

chao1_plot
ggsave('chao1_16s.png', chao1_plot, width = 5, height = 2.5, dpi = 768)

# ---------------- Shannon ----------------
alpha_bac_rare_Shannon <- as.data.frame(alpha_bac_rare[ , "Shannon"])
colnames(alpha_bac_rare_Shannon) <- "Shannon"
alpha_bac_rare_Shannon_final <- cbind(designFinal, alpha_bac_rare_Shannon) %>% as.data.frame()

alpha_bac_rare_Shannon_final$Substrate <- ordered(alpha_bac_rare_Shannon_final$Substrate, levels = c("Chitin","Mineral"))
alpha_bac_rare_Shannon_final$Inoculum <- ordered(alpha_bac_rare_Shannon_final$Inoculum, levels = c("1","3","5","7","9"))
alpha_bac_rare_Shannon_final$Cultivar <- ordered(alpha_bac_rare_Shannon_final$Cultivar, levels = c("Babushka","RGT"))

x.axis.label <- "Microbiome"
y.axis.label <- "Shannon"
alpha_bac_rare_Shannon_final$Inoculum <- as.numeric(levels(alpha_bac_rare_Shannon_final$Inoculum))[alpha_bac_rare_Shannon_final$Inoculum]

shannon_plot <- ggplot(alpha_bac_rare_Shannon_final, aes(Inoculum, Shannon, color = Cultivar, shape = Cultivar)) +
  geom_smooth(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "Babushka"), color = "#ff791f", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "RGT"),      color = "#000000", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "Babushka"), color = "#000000", method = "lm", se = FALSE, linewidth = 1) +
  geom_smooth(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "RGT"),      color = "#FFFFFF", method = "lm", se = FALSE, linewidth = 1) +
  geom_point(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "Babushka"), shape = 16, size = 4, color = "#ff791f") +
  geom_point(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "RGT"),      shape = 17, size = 4, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "Babushka"), shape = 16, size = 3, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "RGT"),      shape = 17, size = 3, color = "#FFFFFF") +
  labs(x = x.axis.label, y = y.axis.label) +
  facet_grid(~ Substrate) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.justification = "right",
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", face = "bold"),
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    legend.text   = element_text(size = 12, color = "black"),
    legend.title  = element_blank(),
    strip.text.x  = element_text(size = 12, color = "black"),
    strip.text.y  = element_text(size = 12, color = "black"),
    panel.spacing = unit(1.5, "lines")
  ) +
  scale_x_continuous(
    breaks = c(1,3,5,7,9),
    labels = c(expression(10^-1), expression(10^-3), expression(10^-5), expression(10^-7), expression("No"))
  ) +
  scale_shape_manual(values = c("Babushka" = 16, "RGT" = 17))

shannon_plot
ggsave('shannon_16s.png', shannon_plot, width = 5, height = 2.5, dpi = 768)

# ---------------- Simpson ----------------
alpha_bac_rare_Simpson <- as.data.frame(alpha_bac_rare[ , "Simpson"])
colnames(alpha_bac_rare_Simpson) <- "Simpson"

alpha_bac_rare_Simpson_final <- cbind(designFinal, alpha_bac_rare_Simpson) %>% as.data.frame()

alpha_bac_rare_Simpson_final$Substrate <- ordered(alpha_bac_rare_Simpson_final$Substrate, levels = c("Chitin","Mineral"))
alpha_bac_rare_Simpson_final$Inoculum <- ordered(alpha_bac_rare_Simpson_final$Inoculum, levels = c("1","3","5","7","9"))
alpha_bac_rare_Simpson_final$Cultivar <- ordered(alpha_bac_rare_Simpson_final$Cultivar, levels = c("Babushka","RGT"))

x.axis.label <- "Microbiome"
y.axis.label <- "Simpson"
alpha_bac_rare_Simpson_final$Inoculum <- as.numeric(levels(alpha_bac_rare_Simpson_final$Inoculum))[alpha_bac_rare_Simpson_final$Inoculum]

simpson_plot <- ggplot(alpha_bac_rare_Simpson_final, aes(Inoculum, Simpson, color = Cultivar, shape = Cultivar)) +
  geom_smooth(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "Babushka"), color = "#ff791f", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "RGT"),      color = "#000000", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "Babushka"), color = "#000000", method = "lm", se = FALSE, linewidth = 1) +
  geom_smooth(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "RGT"),      color = "#FFFFFF", method = "lm", se = FALSE, linewidth = 1) +
  geom_point(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "Babushka"), shape = 16, size = 4, color = "#ff791f") +
  geom_point(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "RGT"),      shape = 17, size = 4, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "Babushka"), shape = 16, size = 3, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "RGT"),      shape = 17, size = 3, color = "#FFFFFF") +
  labs(x = x.axis.label, y = y.axis.label) +
  facet_grid(~ Substrate) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.justification = "right",
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", face = "bold"),
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    legend.text   = element_text(size = 12, color = "black"),
    legend.title  = element_blank(),
    strip.text.x  = element_text(size = 12, color = "black"),
    strip.text.y  = element_text(size = 12, color = "black"),
    panel.spacing = unit(1.5, "lines")
  ) +
  scale_x_continuous(
    breaks = c(1,3,5,7,9),
    labels = c(expression(10^-1), expression(10^-3), expression(10^-5), expression(10^-7), expression("No"))
  ) +
  scale_shape_manual(values = c("Babushka" = 16, "RGT" = 17))

simpson_plot
ggsave('simpson_16s.png', simpson_plot, width = 5, height = 2.5, dpi = 768)

# ---------------------------- #
# 11) Simple linear models (alpha ~ Inoculum) within cultivar × substrate
# ---------------------------- #
data_alphadiv     <- cbind(designFinal, alpha_bac_rare_Observed, alpha_bac_rare_Chao1, alpha_bac_rare_Shannon, alpha_bac_rare_Simpson)
data_alphadiv_df  <- as.data.frame(data_alphadiv)

# Simpson
LM_Simpson_bab_chi <- lm(Simpson ~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Chitin")
LM_Simpson_bab_min <- lm(Simpson ~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Mineral")
LM_Simpson_rgt_chi <- lm(Simpson ~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT"      & Substrate == "Chitin")
LM_Simpson_rgt_min <- lm(Simpson ~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT"      & Substrate == "Mineral")
summary(LM_Simpson_bab_chi); summary(LM_Simpson_bab_min); summary(LM_Simpson_rgt_chi); summary(LM_Simpson_rgt_min)

# Chao1
LM_Chao1_bab_chi <- lm(Chao1 ~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Chitin")
LM_Chao1_bab_min <- lm(Chao1 ~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Mineral")
LM_Chao1_rgt_chi <- lm(Chao1 ~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT"      & Substrate == "Chitin")
LM_Chao1_rgt_min <- lm(Chao1 ~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT"      & Substrate == "Mineral")
summary(LM_Chao1_bab_chi); summary(LM_Chao1_bab_min); summary(LM_Chao1_rgt_chi); summary(LM_Chao1_rgt_min)

# Shannon
LM_Shannon_bab_chi <- lm(Shannon ~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Chitin")
LM_Shannon_bab_min <- lm(Shannon ~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Mineral")
LM_Shannon_rgt_chi <- lm(Shannon ~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT"      & Substrate == "Chitin")
LM_Shannon_rgt_min <- lm(Shannon ~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT"      & Substrate == "Mineral")
summary(LM_Shannon_bab_chi); summary(LM_Shannon_bab_min); summary(LM_Shannon_rgt_chi); summary(LM_Shannon_rgt_min)

# ---------------------------- #
# 12) Figure 6 — Beta diversity (Bray–Curtis) PCoA + PERMANOVA (adonis2)
# ---------------------------- #
# Ordinate with Bray–Curtis on genus-level rarefied 16S data
bac_beta_div_bray_rare <- ordinate(bac_rare_genus, method = "PCoA", distance = "bray")

# Bray–Curtis distance matrix for PERMANOVA (10,000 perms; by = "terms")
BC_bac_clean_rare <- phyloseq::distance(bac_rare_genus, method = "bray")

# Load external metadata for adonis2 (must match sample order / rownames)
# Ensure 'ino' is character for interaction terms
design <- read.delim("metadata.tsv", sep = "\t", header = TRUE, row.names = 1)
design$ino <- as.character(design$ino)

Adonis_bray_combined <- adonis2(BC_bac_clean_rare ~ substrate * ino * plant, data = design,
                                permutations = 10000, by = "terms")

# PCoA plot colored by inoculum, shaped by plant_sub (requires 'plant_sub' in sample_data)
pcoa_plot <- plot_ordination(bac_rare_genus, bac_beta_div_bray_rare, axes = c(1, 2),
                             shape = "plant_sub", color = "ino") +
  geom_point(size = 6, alpha = 5) +
  scale_colour_manual(name = "Inoculum", values = colorblind_Palette) +
  scale_shape_manual(values = c(1, 2, 16, 17)) +
  theme_bw() +
  theme(
    axis.title.y   = element_text(color = "Black", size = 23),
    axis.title.x   = element_text(color = "Black", size = 23),
    legend.key.size= unit(1.5, "line"),
    legend.text    = element_text(size = 18),
    axis.text.y    = element_text(size = 23),
    axis.text.x    = element_text(size = 23),
    legend.position= "right",
    axis.line      = element_line(linewidth = 1),
    axis.ticks     = element_line(linewidth = 1),
    legend.title   = element_text(size = 18, face = "bold")
  )

pcoa_plot
ggsave('PCoA_bac.png', pcoa_plot, width = 7, height = 5, dpi = 768)