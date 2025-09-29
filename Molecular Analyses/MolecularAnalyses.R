# Load libraries
library("DESeq2");packageVersion("DESeq2") 
library ("ape");packageVersion("ape") 
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library (tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr") 
library ("vegan"); packageVersion("vegan")
library("MicrobiotaProcess")
library("microbiome")
library(Polychrome)
library(tidyr)
library(pheatmap)

palette_named <- readRDS("genus_palette.rds")
# OTU tables
otu_its<- read.csv(file = "its_OTU.csv",row.names = 1, sep="\t",header = TRUE)
otu_pro<- read.csv(file = "protist_OTU.csv",row.names = 1, sep="\t",header = TRUE)
otu_16s<- read.csv(file = "16s_OTU.tsv",row.names = 1, sep="\t",header = TRUE)

# Taxonomy tables
taxonomy_its <- read.table(file = "its_taxonomy.tsv", sep = "\t", header = T ,row.names = 1)
taxonomy_pro <- read.table(file ="protist_taxonomy.tsv", sep = "\t", header = T ,row.names = 1)
taxonomy_16s <- read.table(file = "16s_taxonomy.tsv", sep = "\t", header = T ,row.names = 1)

# Clean taxonomy table
# ITS
tax_its <- taxonomy_its %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")
tax_its.clean <- data.frame(row.names = row.names(tax_its),
                            Kingdom = str_replace(tax_its[,1], "d__",""),
                            Phylum = str_replace(tax_its[,2], "p__",""),
                            Class = str_replace(tax_its[,3], "c__",""),
                            Order = str_replace(tax_its[,4], "o__",""),
                            Family = str_replace(tax_its[,5], "f__",""),
                            Genus = str_replace(tax_its[,6], "g__",""),
                            Species = str_replace(tax_its[,7], "s__",""),
                            stringsAsFactors = FALSE)

# Protist
tax_pro <- taxonomy_pro %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")
tax_pro.clean <- data.frame(row.names = row.names(tax_pro),
                            Kingdom = str_replace(tax_pro[,1], "d__",""),
                            Phylum = str_replace(tax_pro[,2], "p__",""),
                            Class = str_replace(tax_pro[,3], "c__",""),
                            Order = str_replace(tax_pro[,4], "o__",""),
                            Family = str_replace(tax_pro[,5], "f__",""),
                            Genus = str_replace(tax_pro[,6], "g__",""),
                            Species = str_replace(tax_pro[,7], "s__",""),
                            stringsAsFactors = FALSE)

# 16S
tax_16s <- taxonomy_16s %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")
tax_16s.clean <- data.frame(row.names = row.names(tax_16s),
                            Kingdom = str_replace(tax_16s[,1], "d__",""),
                            Phylum = str_replace(tax_16s[,2], "p__",""),
                            Class = str_replace(tax_16s[,3], "c__",""),
                            Order = str_replace(tax_16s[,4], "o__",""),
                            Family = str_replace(tax_16s[,5], "f__",""),
                            Genus = str_replace(tax_16s[,6], "g__",""),
                            Species = str_replace(tax_16s[,7], "s__",""),
                            stringsAsFactors = FALSE)

# Metadata files
mapfile_its=import_qiime_sample_data("its_metadata.tsv")
mapfile_pro=import_qiime_sample_data("protist_metadata.tsv")
mapfile_16s=import_qiime_sample_data("16s_metadata.tsv")

# Create to phyloseq object
OTU_ITS = otu_table(as.matrix(otu_its), taxa_are_rows = TRUE)
OTU_PRO = otu_table(as.matrix(otu_pro), taxa_are_rows = TRUE)
OTU_16S = otu_table(as.matrix(otu_16s), taxa_are_rows = TRUE)

TAX_ITS = phyloseq::tax_table(as.matrix(tax_its.clean))
TAX_PRO= phyloseq::tax_table(as.matrix(tax_pro.clean))
TAX_16S = phyloseq::tax_table(as.matrix(tax_16s.clean))

META_ITS <- sample_data(mapfile_its)
META_PRO <- sample_data(mapfile_pro)
META_16S <- sample_data(mapfile_16s)

physeq_its <- phyloseq(OTU_ITS, TAX_ITS, META_ITS)
physeq_pro <- phyloseq(OTU_PRO, TAX_PRO, META_PRO)
physeq_16s <- phyloseq(OTU_16S, TAX_16S, META_16S)

# Clean all three datasets
# ITS
clean_its=subset_taxa(physeq_its, Phylum!= "NA")
clean_its=subset_taxa(clean_its, Genus!= "NA")
clean_its=subset_taxa(clean_its, Genus!= "uncultured")
clean_its=subset_taxa(clean_its, Phylum!= "Fungi_phy_Incertae_sedis")
clean_its=subset_taxa(clean_its, Genus!= "Fungi_gen_Incertae_sedis")
clean_its=subset_taxa(clean_its, Genus!= "")
clean_its=subset_taxa(clean_its, Phylum!= "Streptophyta")
clean_its=subset_taxa(clean_its, Phylum!= "Chlorophyta")
clean_its=subset_taxa(clean_its, Genus!= "Sporidiobolaceae_gen_Incertae_sedis")
clean_its=subset_taxa(clean_its, Genus!= "Lentitheciaceae_gen_Incertae_sedis")
clean_its=subset_taxa(clean_its, Genus!= "GS11_gen_Incertae_sedis")
clean_its=subset_taxa(clean_its, Genus!= "Malasseziaceae_gen_Incertae_sedis")
clean_its=subset_taxa(clean_its, Genus!= "Sebacinales_gen_Incertae_sedis")
clean_its=subset_taxa(clean_its, Genus!= "Candida")
clean_its=subset_taxa(clean_its, (Kingdom=="Fungi") | is.na(Kingdom))
final_clean_its = filter_taxa(clean_its, function(x) sum(x > 1) > (0.005*length(x)), TRUE)
final_clean_its

# Clean protist
clean_pro=subset_taxa(physeq_pro, Phylum!= "NA")
clean_pro=subset_taxa(clean_pro, Genus!= "NA")
clean_pro=subset_taxa(clean_pro, Phylum=="Cercozoa")
clean_pro=subset_taxa(clean_pro, Genus!= "uncultured")
final_clean_pro = filter_taxa(clean_pro, function(x) sum(x > 1) > (0.005*length(x)), TRUE)

# Clean 16S
clean_16s=subset_taxa(physeq_16s, (Order!="Chloroplast") | is.na(Order))
clean_16s=subset_taxa(clean_16s, (Family!="Mitochondria") | is.na(Family))
clean_16s=subset_taxa(clean_16s, Genus!= "NA")
clean_16s=subset_taxa(clean_16s, Phylum!= "NA")
clean_16s=subset_taxa(clean_16s, Genus!= "uncultured")
final_clean_16s = filter_taxa(clean_16s, function(x) sum(x > 1) > (0.005*length(x)), TRUE)

# 16s data
# Rarefy
bac_rare <- rarefy_even_depth(final_clean_16s, rngseed=TRUE)

# Agglomerate
bac_rare_phylum <- tax_glom(bac_rare, taxrank= "Phylum", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
bac_rare_genus <- tax_glom(bac_rare, taxrank= "Genus", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

# Figure 3
# Relative abundance (per treatment)
# Calculate relative abundance 
bac_rel = transform_sample_counts(bac_rare, function(x) x/sum(x)*100)
# agglomerate taxa
bac_rel_genus <- tax_glom(bac_rel, taxrank = "Genus")
bac_rel_genus_melt <- psmelt(bac_rel_genus)
# Summarize mean abundance per treatment and genus
bac_rel_genus_melt_sum <- bac_rel_genus_melt %>%
  group_by(substrate, plant, ino, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = 'drop') %>%
  mutate(Treatment = paste(plant, ino, sep = "_"))

# Identify top 10 genera per treatment group
top10_genera <- bac_rel_genus_melt_sum %>%
  group_by(substrate, plant, ino) %>%
  slice_max(order_by = Abundance, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Top10 = TRUE)

# Join top 10 info back to full data and label non-top10 as "Others"
melt_top10 <- bac_rel_genus_melt_sum %>%
  left_join(top10_genera %>%
              select(substrate, plant, ino, Genus, Top10),
            by = c("substrate", "plant", "ino", "Genus")) %>%
  mutate(Genus = ifelse(is.na(Top10), "Others", as.character(Genus)))

# Collapse "Others" and summarize again
melt_top10 <- melt_top10 %>%
  group_by(substrate, plant, ino, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(Treatment = paste(plant, ino, sep = "_"))

#  Define consistent global Genus levels (Others last)
all_genera <- unique(melt_top10$Genus) %>% as.character()
genus_levels <- c(sort(setdiff(all_genera, "Others")), "Others")
melt_top10$Genus <- factor(melt_top10$Genus, levels = genus_levels)

# Plot
barplot_top10=ggplot(melt_top10, aes(x = ino, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_grid(substrate ~ plant, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = palette_named) +
  labs(#title = "Top 10 Genera per Treatment (Others Collapsed)",
    x = "ino",
    y = "Mean Relative Abundance (%)") +
  theme_minimal() +
  theme(
    legend.position = "left",
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 16, color = "black", face = "bold"),
    strip.text = element_text(size = 16, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.ticks.y = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length.y = unit(0.3, "cm")
  ) +
  labs(y = "Relative Abundance (%)")
barplot_top10
ggsave('barplot_top10.png',barplot_top10,width=12,height=10, dpi = 768)

# Figure 4 
library(randomForest)

# Create wide abundance matrix: samples x genera
abundance_wide <- bac_rel_genus_melt %>%
  ungroup()%>%
  select(sample = Sample, Genus, Abundance) %>%
  pivot_wider(names_from = Genus, values_from = Abundance, values_fill = 0)
colnames(abundance_wide) <- make.names(colnames(abundance_wide))

# Combine with treatment labels
sample_treatment <- bac_rel_genus_melt %>%
  ungroup() %>%
  mutate(Treatment = paste(substrate, plant, ino, sep = "_")) %>%
  select(sample = Sample, Treatment) %>%
  distinct()
rf_data <- left_join(sample_treatment, abundance_wide, by = "sample")

# Random Forest model
set.seed(123)
rf_model <- randomForest(as.factor(Treatment) ~ ., ntree=50000,data = rf_data %>% 
                           select(-sample), importance = TRUE)

# Extract important genera
imp_genera <- importance(rf_model)
imp_genera_df <- data.frame(Genus = rownames(imp_genera), MeanDecreaseGini = imp_genera[, "MeanDecreaseGini"], MeanDecreaseAccuracy=imp_genera[,"MeanDecreaseAccuracy" ])

top_15_important <- imp_genera_df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  slice_head(n = 15)

# Plot
# Subset abundance matrix to top important genera only
top_genera_names <- top_15_important$Genus

# 1. Prepare abundance data by averaging across treatments
abundance_top <- rf_data %>%
  select(sample, Treatment, all_of(top_genera_names))

abundance_avg <- abundance_top %>%
  group_by(Treatment) %>%
  summarise_at(vars(all_of(top_genera_names)), mean) %>%
  ungroup()
# Convert to matrix and set rownames as Treatment
abundance_avg_mat <- abundance_avg %>%
  column_to_rownames(var = "Treatment") %>%
  as.matrix()

# 2. Create annotation for columns (Treatment info)
annotation_df <- data.frame(Treatment = rownames(abundance_avg_mat)) %>%
  mutate(
    substrate = sapply(strsplit(as.character(Treatment), "_"), `[`, 1),
    plant     = sapply(strsplit(as.character(Treatment), "_"), `[`, 2),
    ino       = sapply(strsplit(as.character(Treatment), "_"), `[`, 3)
  ) %>%
  select(-Treatment) %>%
  mutate(
    substrate = factor(substrate),
    plant     = factor(plant),
    ino       = factor(ino)
  )
rownames(annotation_df) <- rownames(abundance_avg_mat)

# 3. Scale abundance matrix by genus (rows)
abundance_scaled <- t(scale(t(abundance_avg_mat)))  # scale by row (genus)

# 4. Create color mappings for annotation bars
substrate_levels <- levels(annotation_df$substrate)
plant_levels     <- levels(annotation_df$plant)
ino_levels       <- levels(annotation_df$ino)

substrate_colors <- c(
  Chitin = "#b892ff",     
  Mineral = "#4cc9f0")    

plant_colors <- c(
  Babushka = "#80ed99",     
  RGT = "#ffafcc")    

ino_colors <- c(
  "10-1" = "#ff686b",   
  "10-3" = "#ff9b85",
  "10-5"="#ffd97d",
  "10-7"="#aaf683",
  "no. ino"="#60d394") 

annotation_colors <- list(
  substrate = substrate_colors,
  plant     = plant_colors,
  ino       = ino_colors
)

annotation_df <- annotation_df %>%
  select(ino, plant, substrate)

# 5. Plot heatmap (transposed so genera are rows and treatments are columns)
heatmap_transposed <- pheatmap(
  t(abundance_scaled),              # transpose: genera are rows
  annotation_col = annotation_df,   # annotate columns (treatments)
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  # display_number=TRUE,
  # display_numbers = round(t(abundance_avg_mat), 2),
  annotation_colors = annotation_colors
)

ggsave('heatmap_zscore_unclustered.png',heatmap_transposed,width =6, height = 3, dpi = 600)


# Figure 5
# Estimate alpha diversity
alpha_bac_rare<-  estimate_richness(bac_rare_genus, measures = c("Observed", "Shannon", "Chao1","Simpson"))

# Design file
design <- bac_rare_genus@sam_data

# substrate
design_substrate=as.data.frame(design[,6])
rownames(design_substrate) =rownames(design)
colnames(design_substrate) <- c("Substrate")

# cultivar
design_cultivar = as.data.frame(design[, 4])
rownames(design_cultivar) <- rownames(design)
colnames(design_cultivar) <- c("Cultivar")

# inoculum
design_inoculum <- as.data.frame(design[, 5])
rownames(design_inoculum) <- rownames(design)
colnames(design_inoculum) <- c("Inoculum")
design_inoculum 

designFinal <- cbind(design_cultivar, design_inoculum, design_substrate)

# Observed
alpha_bac_rare_Observed <- as.data.frame(alpha_bac_rare[ ,1])
rownames(alpha_bac_rare_Observed) <- rownames(alpha_bac_rare)
colnames(alpha_bac_rare_Observed) <- c("Observed")

# Combine
alpha_bac_rare_Observed_final <- cbind(designFinal, alpha_bac_rare)
alpha_bac_rare_Observed_final <- as.data.frame(alpha_bac_rare_Observed_final)

# Reorder
alpha_bac_rare_Observed_final$Substrate <- ordered(alpha_bac_rare_Observed_final$Substrate, levels=c("Chitin","Mineral")) 
alpha_bac_rare_Observed_final$Inoculum <- ordered(alpha_bac_rare_Observed_final$Inoculum, levels=c("1", "3","5","7","9"))  
alpha_bac_rare_Observed_final$Cultivar <- ordered(alpha_bac_rare_Observed_final$Cultivar, levels=c("Babushka","RGT")) 

# Rename
cultivar_names <- c(`Babushka`="Babushka",
                    `RGT`="RGT")
substrate_names <- c(
  `Chitin` = "Chitin",
  `Mineral` = "Mineral")
inoculum_names<-c(`10-1`="1",`10-3`="3",`10-5`="5",`10-7`="7",`no. ino`="9")

# Chao1
alpha_bac_rare_Chao1 <- as.data.frame(alpha_bac_rare[ ,2])
rownames(alpha_bac_rare_Chao1) <- rownames(alpha_bac_rare)
colnames(alpha_bac_rare_Chao1) <- c("Chao1")

# Combine
alpha_bac_rare_Chao1_final <- cbind(designFinal, alpha_bac_rare_Chao1)
alpha_bac_rare_Chao1_final <- as.data.frame(alpha_bac_rare_Chao1_final)

# Reorder
alpha_bac_rare_Chao1_final$Substrate <- ordered(alpha_bac_rare_Chao1_final$Substrate, levels=c("Chitin","Mineral")) 
alpha_bac_rare_Chao1_final$Inoculum <- ordered(alpha_bac_rare_Chao1_final$Inoculum, levels=c("1", "3","5","7","9"))
alpha_bac_rare_Chao1_final$Cultivar <- ordered(alpha_bac_rare_Chao1_final$Cultivar, levels=c("Babushka","RGT")) 


# Plot axis
x.axis.label="Microbiome"
y.axis.label="Chao1"
alpha_bac_rare_Chao1_final$Inoculum=as.numeric(levels(alpha_bac_rare_Chao1_final$Inoculum))[alpha_bac_rare_Chao1_final$Inoculum]

# Plot chao1
chao1_plot=ggplot(alpha_bac_rare_Chao1_final,aes(Inoculum,Chao1,color=Cultivar, shape=Cultivar))+
  geom_smooth(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "Babushka"),color = "#ff791f", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "RGT"), color = "#000000", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "Babushka"),color = "#000000", method = "lm", se = FALSE, linewidth = 1) +
  geom_smooth(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "RGT"), color = "#FFFFFF", method = "lm", se = FALSE, linewidth = 1) +
  geom_point(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "Babushka"), shape = 16, size = 4, color = "#ff791f") +
  geom_point(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "RGT"), shape = 17, size = 4, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "Babushka"), shape = 16, size = 3, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Chao1_final, Cultivar == "RGT"), shape = 17, size = 3, color = "#FFFFFF") +
  labs(x = x.axis.label, y = y.axis.label) +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  facet_grid(~ Substrate) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.justification = "right",
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 12, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 12, color = "black"),
        strip.text.y = element_text(size = 12, color = "black"),
        panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks = c(1, 3,5,7,9),
                     labels = c(expression(10^-1),
                                expression(10^-3),expression(10^-5),expression(10^-7),expression("No"))) +
  scale_shape_manual(values = c("Babushka" = 16, "RGT" = 17))+
  theme(axis.title.y=element_text(size=12,color="black"))+
  theme(legend.text=element_text(size=10,color="black"))+
  theme(legend.title=element_text(size=10,color="black"))

# Display and save the plot
chao1_plot
ggsave('chao1_16s.png',chao1_plot,width=5,height=2.5, dpi = 768)


# Shannon
alpha_bac_rare_Shannon <- as.data.frame(alpha_bac_rare[ ,4])
rownames(alpha_bac_rare_Shannon) <- rownames(alpha_bac_rare)
colnames(alpha_bac_rare_Shannon) <- c("Shannon")

# Combine
alpha_bac_rare_Shannon_final <- cbind(designFinal, alpha_bac_rare_Shannon)
alpha_bac_rare_Shannon_final <- as.data.frame(alpha_bac_rare_Shannon_final)

# Reorder
alpha_bac_rare_Shannon_final$Substrate <- ordered(alpha_bac_rare_Shannon_final$Substrate, levels=c("Chitin","Mineral"))
alpha_bac_rare_Shannon_final$Inoculum <- ordered(alpha_bac_rare_Shannon_final$Inoculum, levels=c("1", "3","5","7","9"))
alpha_bac_rare_Shannon_final$Cultivar <- ordered(alpha_bac_rare_Shannon_final$Cultivar, levels=c("Babushka","RGT"))

# Plot axis
x.axis.label="Microbiome" 
y.axis.label="Shannon"
alpha_bac_rare_Shannon_final$Inoculum=as.numeric(levels(alpha_bac_rare_Shannon_final$Inoculum))[alpha_bac_rare_Shannon_final$Inoculum]

shannon_plot=ggplot(alpha_bac_rare_Shannon_final,aes(Inoculum,Shannon,color=Cultivar, shape=Cultivar))+
  geom_smooth(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "Babushka"),color = "#ff791f", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "RGT"), color = "#000000", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "Babushka"),color = "#000000", method = "lm", se = FALSE, linewidth = 1) +
  geom_smooth(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "RGT"), color = "#FFFFFF", method = "lm", se = FALSE, linewidth = 1) +
  geom_point(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "Babushka"), shape = 16, size = 4, color = "#ff791f") +
  geom_point(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "RGT"), shape = 17, size = 4, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "Babushka"), shape = 16, size = 3, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Shannon_final, Cultivar == "RGT"), shape = 17, size = 3, color = "#FFFFFF") +
  labs(x = x.axis.label, y = y.axis.label) +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  facet_grid(~ Substrate) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.justification = "right",
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 12, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 12, color = "black"),
        strip.text.y = element_text(size = 12, color = "black"),
        panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks = c(1, 3,5,7,9),
                     labels = c(expression(10^-1),
                                expression(10^-3),expression(10^-5),expression(10^-7),expression("No"))) +
  scale_shape_manual(values = c("Babushka" = 16, "RGT" = 17))+
  theme(axis.title.y=element_text(size=12,color="black"))+
  theme(legend.text=element_text(size=10,color="black"))+
  theme(legend.title=element_text(size=10,color="black"))

# Display and save the plot
shannon_plot
ggsave('shannon_16s.png',shannon_plot,width=5,height=2.5, dpi = 768)


# Simpson
alpha_bac_rare_Simpson <- as.data.frame(alpha_bac_rare[ ,5])
rownames(alpha_bac_rare_Simpson) <- rownames(alpha_bac_rare)
colnames(alpha_bac_rare_Simpson) <- c("Simpson")

#Combine the dataset sample description and Simpson Genera
alpha_bac_rare_Simpson_final <- cbind(designFinal, alpha_bac_rare_Simpson)
alpha_bac_rare_Simpson_final <- as.data.frame(alpha_bac_rare_Simpson_final)

#Order the levels according to a defined order
alpha_bac_rare_Simpson_final$Substrate <- ordered(alpha_bac_rare_Simpson_final$Substrate, levels=c("Chitin","Mineral"))
alpha_bac_rare_Simpson_final$Inoculum <- ordered(alpha_bac_rare_Simpson_final$Inoculum, levels=c("1", "3","5","7","9"))
alpha_bac_rare_Simpson_final$Cultivar <- ordered(alpha_bac_rare_Simpson_final$Cultivar, levels=c("Babushka","RGT"))

x.axis.label="Microbiome" #<-- define the label of x axis
y.axis.label="Simpson" #<-- define the label of x axis
alpha_bac_rare_Simpson_final$Inoculum=as.numeric(levels(alpha_bac_rare_Simpson_final$Inoculum))[alpha_bac_rare_Simpson_final$Inoculum]

simpson_plot=ggplot(alpha_bac_rare_Simpson_final,aes(Inoculum,Simpson,color=Cultivar, shape=Cultivar))+
  geom_smooth(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "Babushka"),color = "#ff791f", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "RGT"), color = "#000000", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "Babushka"),color = "#000000", method = "lm", se = FALSE, linewidth = 1) +
  geom_smooth(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "RGT"), color = "#FFFFFF", method = "lm", se = FALSE, linewidth = 1) +
  geom_point(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "Babushka"), shape = 16, size = 4, color = "#ff791f") +
  geom_point(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "RGT"), shape = 17, size = 4, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "Babushka"), shape = 16, size = 3, color = "#000000") +
  geom_point(data = subset(alpha_bac_rare_Simpson_final, Cultivar == "RGT"), shape = 17, size = 3, color = "#FFFFFF") +
  labs(x = x.axis.label, y = y.axis.label) +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  facet_grid(~ Substrate) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.justification = "right",
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 12, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 12, color = "black"),
        strip.text.y = element_text(size = 12, color = "black"),
        panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks = c(1, 3,5,7,9),
                     labels = c(expression(10^-1),
                                expression(10^-3),expression(10^-5),expression(10^-7),expression("No"))) +
  scale_shape_manual(values = c("Babushka" = 16, "RGT" = 17))+
  theme(axis.title.y=element_text(size=12,color="black"))+
  theme(legend.text=element_text(size=10,color="black"))+
  theme(legend.title=element_text(size=10,color="black"))

# Display the plot
simpson_plot
ggsave('simpson_16s.png',simpson_plot,width=5,height=2.5, dpi = 768)

# Linera models
data_alphadiv=cbind(designFinal,alpha_bac_rare_Observed, alpha_bac_rare_Chao1, alpha_bac_rare_Shannon, alpha_bac_rare_Simpson)
data_alphadiv_df=as.data.frame(data_alphadiv)

# Simpson
LM_Simpson_bab_chi <- lm(Simpson~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Chitin")
LM_Simpson_bab_min <- lm(Simpson~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Mineral")
LM_Simpson_rgt_chi <- lm(Simpson~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT" & Substrate == "Chitin")
LM_Simpson_rgt_min <- lm(Simpson~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT" & Substrate == "Mineral")
summary(LM_Simpson_bab_chi)
summary(LM_Simpson_bab_min)
summary(LM_Simpson_rgt_chi)
summary(LM_Simpson_rgt_min)

# Chao1
LM_Chao1_bab_chi <- lm(Chao1~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Chitin")
LM_Chao1_bab_min <- lm(Chao1~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Mineral")
LM_Chao1_rgt_chi <- lm(Chao1~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT" & Substrate == "Chitin")
LM_Chao1_rgt_min <- lm(Chao1~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT" & Substrate == "Mineral")
summary(LM_Chao1_bab_chi)
summary(LM_Chao1_bab_min)
summary(LM_Chao1_rgt_chi)
summary(LM_Chao1_rgt_min)

# Shannon
LM_Shannon_bab_chi <- lm(Shannon~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Chitin")
LM_Shannon_bab_min <- lm(Shannon~ Inoculum, data = data_alphadiv, subset = Cultivar == "Babushka" & Substrate == "Mineral")
LM_Shannon_rgt_chi <- lm(Shannon~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT" & Substrate == "Chitin")
LM_Shannon_rgt_min <- lm(Shannon~ Inoculum, data = data_alphadiv, subset = Cultivar == "RGT" & Substrate == "Mineral")
summary(LM_Shannon_bab_chi)
summary(LM_Shannon_bab_min)
summary(LM_Shannon_rgt_chi)
summary(LM_Shannon_rgt_min)

# Figure 6
# Beta diversity- PCoA (Bacteria)
bac_beta_div_bray_rare <- ordinate(bac_rare_genus, "PCoA", "bray")

# adonis2
BC_bac_clean_rare <- phyloseq::distance(bac_rare_genus, "bray")
design <- read.delim("../combined/input/input_16S/metadata.tsv", sep = "\t", header=TRUE, row.names=1)
design$ino=as.character(design$ino)
Adonis_bray_combined <- adonis2(BC_bac_clean_rare ~ substrate*ino*plant , data= design, permutations = 10000, by="terms")

# Plot
pcoa_plot=plot_ordination(bac_rare_genus, bac_beta_div_bray_rare, axes = c(1, 2), shape ="plant_sub", color = "ino")+ 
  geom_point(size = 6, alpha = 5)+
  scale_colour_manual(name = "Inoculum", values = colorblind_Palette)+
  scale_shape_manual(values = c(1,2,16,17)) + 
  theme_bw()+
  theme(axis.title.y = element_text(color="Black", size=23),
        axis.title.x = element_text(color="Black", size=23),
        legend.key.size = unit(1.5,"line"),
        legend.text = element_text(size = 18),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23),
        legend.position="right",
        axis.line = element_line(linewidth=1),
        axis.ticks = element_line(linewidth=1),
        legend.title = element_text(size = 18, face="bold")
  )
pcoa_plot
ggsave('PCoA_bac.png',pcoa_plot,width=7,height=5, dpi = 768)
