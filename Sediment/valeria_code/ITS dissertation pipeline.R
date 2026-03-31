############################################################
# ITS dissertation pipeline (ONE SCRIPT) — ASV + Metadata + Taxonomy
# - Reads: asv.table.txt, cleaned_metadata.csv, ITS_tax_best.tsv (optional)
# - Builds BestRank (deepest available taxonomic rank)
# - Produces plots labeled with SampleID, grouped by Location + Zone
# - Exports tables + figures to exports/tables and exports/figures
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggrepel)
})

# =========================
# 0) SETUP (YOUR PATH)
# =========================
project_dir <- "C:/Users/valer/OneDrive/Desktop/Amplicon_P.R_Sediments/ITS_R/ITS_ampl_sed"
setwd(project_dir)

out_fig <- file.path(project_dir, "exports", "figures")
out_tab <- file.path(project_dir, "exports", "tables")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

# Inputs
asv_path  <- file.path(project_dir, "asv.table.txt")
meta_path <- file.path(project_dir, "cleaned_metadata.csv")
tax_path  <- file.path(project_dir, "ITS_tax_best.tsv")  # optional but recommended

# Checks
stopifnot(file.exists(asv_path))
stopifnot(file.exists(meta_path))
if(!file.exists(tax_path)) message("⚠️ ITS_tax_best.tsv not found — taxonomy/BestRank steps will be skipped.")

# =========================
# 1) READ METADATA (EXACT COLUMNS)
# =========================
## =========================
# 1) READ METADATA (CLEAN NAMES)
# =========================
library(janitor)

meta_raw <- read_csv(meta_path, show_col_types = FALSE) %>%
  janitor::clean_names()

# Your columns become:
# zone, sample_id, site, site_abbreviation, location, ...

meta <- meta_raw %>%
  transmute(
    SampleID   = as.character(sample_id),
    Location   = as.character(location),
    Zone       = as.character(zone),
    Site       = as.character(site),
    SiteAbbrev = as.character(site_abbreviation)
  ) %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  mutate(
    Zone = factor(Zone, levels = sort(unique(Zone))),
    Location = factor(Location, levels = sort(unique(Location)))
  )
names(meta_raw)
head(meta)


# =========================
# 2) READ ASV TABLE (your asv.table.txt format)
#   Row 1 = feature IDs (ASV sequences)
#   Rows 2+ = counts, col1 = SampleID
# =========================
asv_raw <- read.delim(asv_path, header = FALSE, sep = "\t", check.names = FALSE)

feature_ids <- as.character(unlist(asv_raw[1, ]))
asv_counts  <- asv_raw[-1, , drop = FALSE]

sample_ids <- as.character(asv_counts[[1]])
otu_mat <- asv_counts[, -1, drop = FALSE] %>% as.matrix()
mode(otu_mat) <- "numeric"

rownames(otu_mat) <- sample_ids
colnames(otu_mat) <- feature_ids

# =========================
# 3) ALIGN SAMPLES WITH METADATA
# =========================
common_samples <- intersect(rownames(otu_mat), meta$SampleID)
if(length(common_samples) == 0) stop("No overlapping SampleID between ASV table and metadata.")

otu_mat <- otu_mat[common_samples, , drop = FALSE]
meta <- meta %>% filter(SampleID %in% common_samples) %>% distinct(SampleID, .keep_all = TRUE)
meta <- meta %>% slice(match(rownames(otu_mat), SampleID))
stopifnot(identical(meta$SampleID, rownames(otu_mat)))

# Order samples by Zone, Location, SampleID for consistent x-axis
sample_order <- meta %>%
  arrange(Zone, Location, SampleID) %>%
  pull(SampleID)

# =========================
# 4) TAXONOMY: BestRank + Depth Summary (optional)
# =========================
tax_best_rank <- NULL
depth_summary <- NULL

if(file.exists(tax_path)){
  tax_best <- read_tsv(tax_path, show_col_types = FALSE) %>%
    rename(FeatureID = `Feature ID`)
  
  tax_best_rank <- tax_best %>%
    mutate(BestRank = case_when(
      !is.na(Species) ~ Species,
      !is.na(Genus)   ~ Genus,
      !is.na(Family)  ~ Family,
      !is.na(Order)   ~ Order,
      !is.na(Class)   ~ Class,
      !is.na(Phylum)  ~ Phylum,
      !is.na(Kingdom) ~ Kingdom,
      TRUE ~ "Unassigned"
    ))
  
  depth_summary <- tax_best_rank %>%
    summarise(
      species      = sum(!is.na(Species)),
      genus        = sum(is.na(Species) & !is.na(Genus)),
      family       = sum(is.na(Genus) & !is.na(Family)),
      order        = sum(is.na(Family) & !is.na(Order)),
      class        = sum(is.na(Order) & !is.na(Class)),
      phylum       = sum(is.na(Class) & !is.na(Phylum)),
      kingdom_only = sum(is.na(Phylum) & !is.na(Kingdom)),
      unassigned   = sum(is.na(Kingdom))
    )
  
  write_csv(depth_summary, file.path(out_tab, "TABLE_Taxonomy_BestRank_DepthSummary.csv"))
}

# =========================
# 5) BUILD LONG TABLE (ASV x sample) + merge metadata + taxonomy
# =========================
asv_long <- as.data.frame(otu_mat) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "FeatureID", values_to = "Abundance") %>%
  left_join(meta, by = "SampleID")

if(!is.null(tax_best_rank)){
  asv_long <- asv_long %>% left_join(tax_best_rank, by = c("FeatureID" = "FeatureID"))
} else {
  asv_long <- asv_long %>% mutate(BestRank = FeatureID)
}

# Ensure BestRank exists even if taxonomy missing/partial
asv_long <- asv_long %>%
  mutate(BestRank = case_when(
    !is.na(BestRank) ~ BestRank,
    TRUE ~ "Unassigned"
  ))

write_csv(asv_long, file.path(out_tab, "TABLE_ASV_long_with_meta_tax.csv"))

# Collapse to BestRank (sum ASVs sharing the same deepest name)
best_rank_long <- asv_long %>%
  group_by(SampleID, Location, Zone, BestRank) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

best_rank_wide <- best_rank_long %>%
  pivot_wider(names_from = BestRank, values_from = Abundance, values_fill = 0)

write_csv(best_rank_wide, file.path(out_tab, "TABLE_BestRank_wide_counts.csv"))

# 
# Drop zero-read samples (prevents NA distances/ordination coords)
keep <- rowSums(otu_mat) > 0
if(any(!keep)){
  message("⚠️ Dropping ", sum(!keep), " samples with 0 total reads: ",
          paste(rownames(otu_mat)[!keep], collapse = ", "))
}
otu_mat <- otu_mat[keep, , drop = FALSE]
meta <- meta %>% filter(SampleID %in% rownames(otu_mat)) %>%
  slice(match(rownames(otu_mat), SampleID))
stopifnot(identical(meta$SampleID, rownames(otu_mat)))

# Recompute rel_mat after filtering
rel_mat <- sweep(otu_mat, 1, rowSums(otu_mat), "/")
rel_mat[is.na(rel_mat)] <- 0


# =========================
# 6) NORMALIZATION (relative abundance) for ordination
# =========================
rel_mat <- sweep(otu_mat, 1, rowSums(otu_mat), "/")
rel_mat[is.na(rel_mat)] <- 0

# =========================
# 7) ALPHA DIVERSITY (from counts) — labeled by SampleID
# =========================
alpha_df <- tibble(
  SampleID = rownames(otu_mat),
  Observed = rowSums(otu_mat > 0),
  Shannon  = vegan::diversity(otu_mat, index = "shannon")
) %>%
  left_join(meta, by = "SampleID") %>%
  mutate(
    SampleID = factor(SampleID, levels = sample_order)
  )

write_csv(alpha_df, file.path(out_tab, "TABLE_AlphaDiversity_fromCounts.csv"))

p_alpha_obs <- ggplot(alpha_df, aes(x = Location, y = Observed, color = Zone)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 2) +
  ggrepel::geom_text_repel(aes(label = as.character(SampleID)), size = 3, max.overlaps = 80) +
  theme_bw() +
  labs(title = "ITS Alpha Diversity (Observed ASVs) by Location", x = "Location", y = "Observed ASVs")

ggsave(file.path(out_fig, "FIG_Alpha_Observed_byLocation_labels.pdf"),
       p_alpha_obs, width = 12, height = 6)

p_alpha_shan <- ggplot(alpha_df, aes(x = Location, y = Shannon, color = Zone)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 2) +
  ggrepel::geom_text_repel(aes(label = as.character(SampleID)), size = 3, max.overlaps = 80) +
  theme_bw() +
  labs(title = "ITS Alpha Diversity (Shannon) by Location", x = "Location", y = "Shannon")

ggsave(file.path(out_fig, "FIG_Alpha_Shannon_byLocation_labels.pdf"),
       p_alpha_shan, width = 12, height = 6)

ggsave(file.path(out_fig, "FIG_Alpha_Observed_byLocation_facetZone_labels.pdf"),
       p_alpha_obs + facet_wrap(~Zone, scales = "free_x"), width = 13, height = 7)

ggsave(file.path(out_fig, "FIG_Alpha_Shannon_byLocation_facetZone_labels.pdf"),
       p_alpha_shan + facet_wrap(~Zone, scales = "free_x"), width = 13, height = 7)


pcoa_df %>% filter(is.na(PCoA1) | is.na(PCoA2) | is.na(Location) | is.na(Zone)) %>% count()
pcoa_df %>% filter(is.na(PCoA1) | is.na(PCoA2) | is.na(Location) | is.na(Zone)) %>% pull(SampleID)

# =========================
# 8) BETA DIVERSITY — PCoA + NMDS (Bray-Curtis) labeled by SampleID
# =========================

bray2 <- as.matrix(vegdist(rel_mat, method = "bray"))
summary(as.vector(bray2[upper.tri(bray2)]))
sum(bray2[upper.tri(bray2)] == 0)
# Hellinger transform
hel_mat <- vegan::decostand(otu_mat, method = "hellinger")

bray_hel <- vegan::vegdist(hel_mat, method = "bray")

# ---- PCoA (recommended final version)
pcoa_hel <- cmdscale(bray_hel, eig = TRUE, k = 2)

pcoa_hel_df <- as.data.frame(pcoa_hel$points) %>%
  setNames(c("PCoA1","PCoA2")) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta, by = "SampleID")

p_pcoa_hel <- ggplot(pcoa_hel_df,
                     aes(PCoA1, PCoA2,
                         color = Location,
                         shape = Zone)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = SampleID),
                           size = 3,
                           max.overlaps = 80) +
  theme_bw() +
  labs(title = "ITS PCoA (Bray-Curtis, Hellinger transformed)",
       x = "PCoA1",
       y = "PCoA2")

ggsave(file.path(out_fig,
                 "FIG_PCoA_Bray_Hellinger_labels.pdf"),
       p_pcoa_hel,
       width = 11,
       height = 7)

# SKIP---- NMDS
set.seed(1)
nmds <- vegan::metaMDS(rel_mat, distance = "bray", k = 2, trymax = 200, autotransform = FALSE)
nmds_df <- as.data.frame(vegan::scores(nmds, display = "sites")) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta, by = "SampleID")

p_nmds <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = Location, shape = Zone)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = 80) +
  theme_bw() +
  labs(title = "ITS NMDS (Bray-Curtis)", x = "NMDS1", y = "NMDS2")

ggsave(file.path(out_fig, "FIG_NMDS_Bray_labels.pdf"), p_nmds, width = 11, height = 7)
ggsave(file.path(out_fig, "FIG_NMDS_Bray_facetZone_labels.pdf"),
       p_nmds + facet_wrap(~Zone), width = 13, height = 7)




# =========================
# 9) TOP TAXA STACKED BAR (BestRank) by Zone + Location
# =========================
topN <- 15

best_rank_rel <- best_rank_long %>%
  group_by(SampleID) %>%
  mutate(RelAbund = Abundance / sum(Abundance)) %>%
  ungroup()

top_taxa <- best_rank_rel %>%
  group_by(BestRank) %>%
  summarise(Total = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = topN) %>%
  pull(BestRank)

best_rank_rel2 <- best_rank_rel %>%
  mutate(BestRank2 = if_else(BestRank %in% top_taxa, BestRank, "Other")) %>%
  group_by(SampleID, Location, Zone, BestRank2) %>%
  summarise(RelAbund = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%
  mutate(SampleID = factor(SampleID, levels = sample_order))

p_bar <- ggplot(best_rank_rel2, aes(x = SampleID, y = RelAbund, fill = BestRank2)) +
  geom_col() +
  facet_grid(Zone ~ Location, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = paste0("ITS Top ", topN, " BestRank taxa (Relative Abundance)"),
    x = "SampleID", y = "Relative abundance", fill = "BestRank"
  )

ggsave(file.path(out_fig, "FIG_TopTaxa_BestRank_StackedBar_byZoneLocation.pdf"),
       p_bar, width = 16, height = 9)

# =========================
# 10) SAVE ORDINATIONS + STATS TABLES (optional but handy)
# =========================
write_csv(pcoa_df, file.path(out_tab, "TABLE_PCoA_Bray_coordinates.csv"))
write_csv(nmds_df, file.path(out_tab, "TABLE_NMDS_Bray_coordinates.csv"))

message("✅ DONE!")
message("Figures: ", out_fig)
message("Tables:  ", out_tab)


