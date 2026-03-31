###############################################
# 16S GENUS PIPELINE — Vieques Sediments (FINAL)
# Replicates kept separate (site_abbrev: BBSS-SF_1, _2, _3)
#
# INPUTS:
#  - asv_table.csv                 (samples x sequences or transposed)
#  - asv_id_map_with_taxonomy.csv  (Sequence -> taxonomy ranks)
#  - Meta data  - Copy of sediment toxicity .csv (metadata with sampleid, site_abbrev, zone, bac_seq)
#
# OUTPUTS:
#  - FIGURES/*.png + *.pdf
#  - TABLES/*.csv
###############################################

### =========================
### 0) SETTINGS
### =========================
base_dir <- "C:/Users/valer/OneDrive/Desktop/Amplicon_P.R_Sediments/16S _R/exports/ALL"
setwd(base_dir)

fn_asv_table <- "asv_table.csv"
fn_tax_map   <- "asv_id_map_with_taxonomy.csv"
fn_meta      <- "Meta data  - Copy of sediment toxicity .csv"

out_fig <- file.path(base_dir, "FIGURES")
out_tab <- file.path(base_dir, "TABLES")
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

### =========================
### 1) PACKAGES
### =========================
pkgs <- c("tidyverse","phyloseq","vegan","indicspecies","janitor","stringr","scales","patchwork")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))
theme_set(theme_bw(base_size = 12))

### =========================
### 2) HELPERS
### =========================
pick_col <- function(df, candidates){
  hit <- candidates[candidates %in% names(df)]
  if(length(hit) == 0) return(NA_character_)
  hit[[1]]
}

is_dna_string <- function(x){
  grepl("^[ACGT]+$", x) & nchar(x) > 80
}

maybe_transpose_counts <- function(m, sampleid_vec = NULL){
  rn <- rownames(m); cn <- colnames(m)
  rn_dna <- mean(is_dna_string(rn), na.rm = TRUE)
  cn_dna <- mean(is_dna_string(cn), na.rm = TRUE)
  
  rn_is_sample <- if(!is.null(sampleid_vec)) mean(rn %in% sampleid_vec) else 0
  cn_is_sample <- if(!is.null(sampleid_vec)) mean(cn %in% sampleid_vec) else 0
  
  if((rn_dna > cn_dna) || (cn_is_sample > rn_is_sample)){
    message("🔁 Transposing count table (sequences appear to be rows).")
    return(t(m))
  } else {
    message("✅ Count table orientation looks correct (samples rows, sequences cols).")
    return(m)
  }
}

map_count_samples_to_meta <- function(count_samples, meta2){
  meta_sampleid <- as.character(meta2$sampleid)
  meta_bacseq   <- as.character(meta2$bac_seq)
  
  norm_fastq <- function(x){
    x %>%
      str_replace("\\.fastq\\.gz$","") %>%
      str_replace("\\.fastq$","") %>%
      str_replace("\\.fq\\.gz$","") %>%
      str_replace("\\.fq$","") %>%
      str_replace("_R[12]_001$","") %>%
      str_replace("_L\\d{3}$","") %>%
      str_replace("_L\\d{3}_R[12]$","") %>%
      str_replace("_L\\d{3}_R[12]_001$","") %>%
      str_replace("_S\\d+$","") %>%
      str_trim()
  }
  
  cs_norm <- norm_fastq(count_samples)
  ms_norm <- norm_fastq(meta_sampleid)
  mb_norm <- norm_fastq(meta_bacseq)
  
  idx1 <- match(count_samples, meta_sampleid)  # exact sampleid
  idx2 <- match(cs_norm, ms_norm)              # normalized sampleid
  idx3 <- match(count_samples, meta_bacseq)    # exact bac_seq
  idx4 <- match(cs_norm, mb_norm)              # normalized bac_seq
  
  mapped_idx <- ifelse(!is.na(idx1), idx1,
                       ifelse(!is.na(idx2), idx2,
                              ifelse(!is.na(idx3), idx3,
                                     ifelse(!is.na(idx4), idx4, NA))))
  
  mapped_sampleid <- ifelse(is.na(mapped_idx), NA, meta_sampleid[mapped_idx])
  
  tibble(
    count_sample = count_samples,
    mapped_sampleid = mapped_sampleid,
    map_method = case_when(
      !is.na(idx1) ~ "exact_sampleid",
      is.na(idx1) & !is.na(idx2) ~ "norm_sampleid",
      is.na(idx1) & is.na(idx2) & !is.na(idx3) ~ "exact_bacseq",
      is.na(idx1) & is.na(idx2) & is.na(idx3) & !is.na(idx4) ~ "norm_bacseq",
      TRUE ~ "unmapped"
    )
  )
}

get_col_ci <- function(df, target){
  hit <- names(df)[tolower(trimws(names(df))) == tolower(target)]
  if(length(hit) == 0) return(NA_character_)
  hit[1]
}

### =========================
### 3) LOAD + STANDARDIZE METADATA (clean_names OK)
### =========================
meta_raw <- read.csv(fn_meta, check.names = FALSE, stringsAsFactors = FALSE)
meta <- janitor::clean_names(meta_raw)

sampleid_candidates <- names(meta)[str_detect(names(meta), "sample") & str_detect(names(meta), "id")]
col_sampleid <- if(length(sampleid_candidates) > 0) sampleid_candidates[1] else NA_character_

col_abbrev <- pick_col(meta, c("abbreviation","abbr","site_abbreviation"))
col_zone   <- pick_col(meta, c("zone"))
col_bacseq <- pick_col(meta, c("bac_seq","bacseq"))

if(is.na(col_sampleid)) stop("Metadata missing SampleID-like column (needs something like sampleid).")
if(is.na(col_abbrev))   stop("Metadata missing Abbreviation (expected abbreviation/site_abbreviation).")
if(is.na(col_zone))     stop("Metadata missing Zone (expected zone).")
if(is.na(col_bacseq))   message("⚠️ bac_seq not found; mapping will use sampleid only.")

meta2 <- meta %>%
  mutate(
    sampleid    = as.character(.data[[col_sampleid]]),
    site_abbrev = as.character(.data[[col_abbrev]]),   # keep replicates separate
    zone        = trimws(as.character(.data[[col_zone]])),
    bac_seq     = if(!is.na(col_bacseq)) as.character(.data[[col_bacseq]]) else sampleid
  ) %>%
  distinct(sampleid, .keep_all = TRUE)

message("Metadata columns used:")
message("  sampleid:    ", col_sampleid)
message("  site_abbrev: ", col_abbrev)
message("  zone:        ", col_zone)
message("  bac_seq:     ", ifelse(is.na(col_bacseq), "(none)", col_bacseq))

### =========================
### 4) LOAD ASV TABLE (DO NOT clean_names!)
### =========================
asv_raw <- read.csv(fn_asv_table, check.names = FALSE, stringsAsFactors = FALSE)
names(asv_raw) <- trimws(names(asv_raw))  # only trim, preserve sequences

sampleid_col_asv <- names(asv_raw)[str_detect(tolower(names(asv_raw)), "sample") &
                                     str_detect(tolower(names(asv_raw)), "id")]
if(length(sampleid_col_asv) == 0) stop("ASV table: could not find a SampleID column.")
sampleid_col_asv <- sampleid_col_asv[1]

count_df <- asv_raw %>%
  mutate(sampleid = as.character(.data[[sampleid_col_asv]])) %>%
  select(-all_of(sampleid_col_asv)) %>%
  relocate(sampleid)

count_mat <- count_df %>%
  column_to_rownames("sampleid") %>%
  as.matrix()
mode(count_mat) <- "numeric"

count_mat <- maybe_transpose_counts(count_mat, sampleid_vec = meta2$sampleid)

# Drop all-zero sequences
count_mat <- count_mat[, colSums(count_mat) > 0, drop = FALSE]
message("ASV count_mat dims (after transpose/clean): ", paste(dim(count_mat), collapse=" x "))
message("Example ASV sequence colnames: ", paste(head(colnames(count_mat)), collapse=" | "))

# Map ASV rownames to metadata sampleid if needed
map_tbl <- map_count_samples_to_meta(rownames(count_mat), meta2)
write.csv(map_tbl, file.path(out_tab, "TABLE_SampleName_Mapping_ASV_to_Metadata.csv"), row.names = FALSE)

mapped <- map_tbl %>% filter(map_method != "unmapped")
if(nrow(mapped) == 0) stop("No ASV samples mapped to metadata sampleid. Check TABLE_SampleName_Mapping_ASV_to_Metadata.csv")

count_mat2 <- count_mat[mapped$count_sample, , drop = FALSE]
rownames(count_mat2) <- mapped$mapped_sampleid

if(any(duplicated(rownames(count_mat2)))) {
  message("⚠️ Duplicate sample IDs after mapping — summing duplicates.")
  count_mat2 <- rowsum(count_mat2, group = rownames(count_mat2))
}
count_mat <- count_mat2

### =========================
### 5) LOAD TAXONOMY MAP (DO NOT clean_names!)
### =========================
taxmap_raw <- read.csv(fn_tax_map, check.names = FALSE, stringsAsFactors = FALSE)
names(taxmap_raw) <- trimws(names(taxmap_raw))  # preserve sequences

seq_col <- get_col_ci(taxmap_raw, "Sequence")
if(is.na(seq_col)) stop("Tax map: could not find Sequence column.")

k_col <- get_col_ci(taxmap_raw, "Kingdom")
p_col <- get_col_ci(taxmap_raw, "Phylum")
c_col <- get_col_ci(taxmap_raw, "Class")
o_col <- get_col_ci(taxmap_raw, "Order")
f_col <- get_col_ci(taxmap_raw, "Family")
g_col <- get_col_ci(taxmap_raw, "Genus")

tax_tbl <- taxmap_raw %>%
  transmute(
    sequence = as.character(.data[[seq_col]]),
    kingdom  = if(!is.na(k_col)) as.character(.data[[k_col]]) else "Unassigned",
    phylum   = if(!is.na(p_col)) as.character(.data[[p_col]]) else "Unassigned",
    class    = if(!is.na(c_col)) as.character(.data[[c_col]]) else "Unassigned",
    order    = if(!is.na(o_col)) as.character(.data[[o_col]]) else "Unassigned",
    family   = if(!is.na(f_col)) as.character(.data[[f_col]]) else "Unassigned",
    genus    = if(!is.na(g_col)) as.character(.data[[g_col]]) else "Unassigned"
  ) %>%
  mutate(across(c(kingdom,phylum,class,order,family,genus),
                ~ifelse(is.na(.x) | .x=="", "Unassigned", .x))) %>%
  filter(!is.na(sequence), sequence != "") %>%
  distinct(sequence, .keep_all = TRUE) %>%
  column_to_rownames("sequence") %>%
  as.matrix()

shared_seqs <- intersect(colnames(count_mat), rownames(tax_tbl))
message("Shared sequences count: ", length(shared_seqs))
if(length(shared_seqs) == 0) stop("No shared sequences between ASV table and taxonomy map.")

count_mat <- count_mat[, shared_seqs, drop = FALSE]
tax_tbl   <- tax_tbl[shared_seqs, , drop = FALSE]

### =========================
### 6) ALIGN METADATA TO COUNTS
### =========================
meta_ps <- meta2 %>%
  filter(sampleid %in% rownames(count_mat)) %>%
  column_to_rownames("sampleid")

count_mat <- count_mat[rownames(meta_ps), , drop = FALSE]

message("✅ Final matched samples: ", nrow(count_mat))
message("✅ Final matched sequences: ", ncol(count_mat))

### =========================
### 7) PHYLOSEQ + GENUS (REL ABUND) + RENAME TAXA TO GENUS
### =========================
ps <- phyloseq(
  otu_table(t(count_mat), taxa_are_rows = TRUE),
  tax_table(tax_tbl),
  sample_data(meta_ps)
)

ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

ps_rel   <- transform_sample_counts(ps, function(x) x / sum(x))
ps_genus <- tax_glom(ps_rel, taxrank = "genus", NArm = FALSE)


# ---- Rename taxa IDs to genus labels (fix OTU/tax_table mismatch) ----
tax_df <- as.data.frame(tax_table(ps_genus))

tax_df$genus  <- as.character(tax_df$genus)
tax_df$family <- as.character(tax_df$family)

tax_df$genus_clean <- ifelse(
  is.na(tax_df$genus) | tax_df$genus=="" | tax_df$genus=="Unassigned",
  paste0("Unassigned_", ifelse(is.na(tax_df$family) | tax_df$family=="", "Family", tax_df$family)),
  tax_df$genus
)

# Make unique genus labels
genus_unique <- make.unique(tax_df$genus_clean)

# Save current OTU IDs in the same order as tax_df rows
old_ids <- rownames(tax_df)

# 1) Rename OTU/taxa IDs in phyloseq
taxa_names(ps_genus) <- genus_unique

# 2) Update taxonomy table rownames to match the new taxa IDs
rownames(tax_df) <- genus_unique

# 3) Reassign cleaned genus column and put back into ps_genus
tax_df$genus <- tax_df$genus_clean
tax_table(ps_genus) <- tax_table(as.matrix(tax_df[, c("kingdom","phylum","class","order","family","genus")]))


all(taxa_names(ps_genus) == rownames(tax_table(ps_genus)))
head(taxa_names(ps_genus))

# Should NOT be DNA strings anymore:
sum(grepl("^[ACGT]+$", taxa_names(ps_genus)) & nchar(taxa_names(ps_genus)) > 80)
taxa_names(ps_genus) <- str_replace_all(taxa_names(ps_genus), "Unassigned_", "Unk_")


# IMPORTANT: rename taxa IDs to genus names so plots/tables show genus (not sequences)
tax_df <- as.data.frame(tax_table(ps_genus))
tax_df$genus  <- as.character(tax_df$genus)
tax_df$family <- as.character(tax_df$family)

tax_df$genus_clean <- ifelse(
  is.na(tax_df$genus) | tax_df$genus=="" | tax_df$genus=="Unassigned",
  paste0("Unassigned_", ifelse(is.na(tax_df$family) | tax_df$family=="", "Family", tax_df$family)),
  tax_df$genus
)

genus_unique <- make.unique(tax_df$genus_clean)
taxa_names(ps_genus) <- genus_unique

tax_df$genus <- tax_df$genus_clean
tax_table(ps_genus) <- tax_table(as.matrix(tax_df[, c("kingdom","phylum","class","order","family","genus")]))

### =========================
### 8) ZONE COLORS (YOUR EXACT HEX)
### =========================
zone_palette <- c(
  "1" = "#E5B800",
  "2" = "#009E73",
  "3" = "#E69F00",
  "4" = "#D55E00"
)
sample_data(ps_genus)$zone <- trimws(as.character(sample_data(ps_genus)$zone))

### =========================
### 9) GENUS HEATMAP BY site_abbrev + ZONE BAR (Top 60 by variance)
### =========================
sd <- data.frame(sample_data(ps_genus)) %>%
  rownames_to_column("sampleid") %>%
  mutate(
    site_abbrev = as.character(site_abbrev),
    zone = as.character(zone)
  )

# Order sites by zone then name
site_zone <- sd %>%
  group_by(site_abbrev) %>%
  summarize(zone_mode = names(sort(table(zone), decreasing = TRUE))[1], .groups="drop") %>%
  mutate(zone_mode = factor(zone_mode, levels = names(zone_palette))) %>%
  arrange(zone_mode, site_abbrev)

site_order <- site_zone$site_abbrev

# Genus x sample matrix (already genus-labeled)
otu_genus <- as.data.frame(otu_table(ps_genus))
otu_genus$Genus <- rownames(otu_genus)

otu_long <- otu_genus %>%
  pivot_longer(-Genus, names_to = "sampleid", values_to = "rel_abund") %>%
  left_join(sd %>% select(sampleid, site_abbrev, zone), by = "sampleid")

# Mean genus rel abundance per replicate-site
site_genus <- otu_long %>%
  group_by(site_abbrev, Genus) %>%
  summarize(mean_rel = mean(rel_abund, na.rm = TRUE), .groups="drop")

site_wide <- site_genus %>%
  pivot_wider(names_from = site_abbrev, values_from = mean_rel, values_fill = 0)

mat_site <- site_wide %>%
  column_to_rownames("Genus") %>%
  as.matrix()

# ---- TOP 60 genera by variance across replicate-sites ----
site_var <- apply(mat_site, 1, var, na.rm = TRUE)
topN <- names(sort(site_var, decreasing = TRUE))[1:min(60, length(site_var))]
mat_prev <- mat_site[topN, intersect(site_order, colnames(mat_site)), drop = FALSE]

heat_df <- as.data.frame(mat_prev) %>%
  rownames_to_column("Genus") %>%
  pivot_longer(-Genus, names_to = "Site", values_to = "MeanRel") %>%
  mutate(
    Site  = factor(Site, levels = colnames(mat_prev)),
    Genus = factor(Genus, levels = rev(rownames(mat_prev))),
    log10p = log10(MeanRel + 1e-6)
  )

p_heat <- ggplot(heat_df, aes(x = Site, y = Genus, fill = log10p)) +
  geom_tile() +
  labs(
    title = "Top 60 Genera by Variance Across Sites (replicates separate)",
    subtitle = "Fill = log10(mean rel abundance + 1e-6)",
    x = "Site (replicate)", y = "Genus", fill = "log10(rel+1e-6)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

zone_bar <- site_zone %>%
  filter(site_abbrev %in% colnames(mat_prev)) %>%
  mutate(site_abbrev = factor(site_abbrev, levels = colnames(mat_prev))) %>%
  ggplot(aes(x = site_abbrev, y = 1, fill = zone_mode)) +
  geom_tile() +
  scale_fill_manual(values = zone_palette, drop = FALSE) +
  theme_void() +
  labs(fill = "Zone")

p_heat_combined <- zone_bar / p_heat + plot_layout(heights = c(0.6, 10))


ggsave(file.path(out_fig, "FIG_16S_GenusHeatmap_Top60Var_bySiteAbbrev_withZoneBar.png"),
       p_heat_combined, width = 14, height = 11, dpi = 300)
ggsave(file.path(out_fig, "FIG_16S_GenusHeatmap_Top60Var_bySiteAbbrev_withZoneBar.pdf"),
       p_heat_combined, width = 14, height = 11)

### =========================
### 10) TOP 10 GENERA PER site_abbrev
### =========================
top10_site <- site_genus %>%
  group_by(site_abbrev) %>%
  slice_max(order_by = mean_rel, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(site_abbrev = factor(site_abbrev, levels = site_order))

write.csv(top10_site, file.path(out_tab, "TABLE_Top10Genera_bySiteAbbrev.csv"), row.names = FALSE)

p_top10 <- top10_site %>%
  ggplot(aes(x = reorder(Genus, mean_rel), y = mean_rel)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ site_abbrev, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Top 10 Genera per Site (replicates separate)",
       x = "Genus", y = "Mean relative abundance") +
  theme(strip.text = element_text(size = 9),
        panel.grid.major.y = element_blank())

ggsave(file.path(out_fig, "FIG_16S_Top10Genera_bySiteAbbrev.png"),
       p_top10, width = 16, height = 10, dpi = 300)
ggsave(file.path(out_fig, "FIG_16S_Top10Genera_bySiteAbbrev.pdf"),
       p_top10, width = 16, height = 10)

### =========================
### 11) INDICATOR GENERA (IndVal) — BY site_abbrev
### =========================
otu_g <- as.data.frame(t(otu_table(ps_genus)))  # samples x genus rel
otu_g <- otu_g[rownames(meta_ps), , drop = FALSE]
otu_g_hel <- vegan::decostand(otu_g, method = "hellinger")

group_site <- as.factor(sample_data(ps_genus)$site_abbrev)
names(group_site) <- sample_names(ps_genus)

set.seed(123)
ind <- indicspecies::multipatt(otu_g_hel, group_site, func = "IndVal.g", control = how(nperm = 999))

ind_tbl <- as.data.frame(ind$sign) %>%
  rownames_to_column("Genus") %>%
  as_tibble()

write.csv(ind_tbl, file.path(out_tab, "TABLE_IndVal_Genus_bySiteAbbrev.csv"), row.names = FALSE)

stat_col <- pick_col(ind_tbl, c("stat","indval","s.value","value"))
if(!is.na(stat_col)) {
  ind_plot_df <- ind_tbl %>%
    arrange(desc(.data[[stat_col]])) %>%
    slice_head(n = 30) %>%
    mutate(Genus = factor(Genus, levels = rev(Genus)))

  p_ind <- ggplot(ind_plot_df, aes(x = Genus, y = .data[[stat_col]])) +
    geom_col() +
    coord_flip() +
    labs(title = "Top Indicator Genera (IndVal) by replicate site_abbrev",
         subtitle = "Top 30 genera ranked by IndVal statistic",
         x = "Genus", y = "IndVal statistic") +
    theme(panel.grid.major.y = element_blank())

  ggsave(file.path(out_fig, "FIG_16S_IndVal_Top30Genera_bySiteAbbrev.png"),
         p_ind, width = 10, height = 8, dpi = 300)
  ggsave(file.path(out_fig, "FIG_16S_IndVal_Top30Genera_bySiteAbbrev.pdf"),
         p_ind, width = 10, height = 8)
}

### =========================
### 12) GENUS ↔ METALS / ORGANIC MATTER (Spearman + FDR)
### =========================
meta_df <- data.frame(sample_data(ps_genus)) %>% rownames_to_column("sampleid")

col_org <- pick_col(meta_df, c("percent_org_total"))
col_cu  <- pick_col(meta_df, c("cu_mg_kg"))
col_ni  <- pick_col(meta_df, c("ni_mg_kg"))
col_zn  <- pick_col(meta_df, c("zn_mg_kg"))
col_co  <- pick_col(meta_df, c("co_mg_kg"))

env_cols <- c(col_org, col_cu, col_ni, col_zn, col_co)
env_cols <- env_cols[!is.na(env_cols)]
if(length(env_cols) == 0) stop("No env columns found (organic/metals). Check metadata column names.")

env <- meta_df %>%
  select(sampleid, all_of(env_cols)) %>%
  mutate(across(-sampleid, as.numeric)) %>%
  column_to_rownames("sampleid")

otu_rel <- as.data.frame(t(otu_table(ps_genus))) # samples x genus rel
otu_rel <- otu_rel[rownames(env), , drop = FALSE]

cor_results <- expand.grid(
  Genus = colnames(otu_rel),
  Variable = colnames(env),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    test = list(suppressWarnings(cor.test(otu_rel[[Genus]], env[[Variable]], method = "spearman", exact = FALSE))),
    rho = unlist(test$estimate),
    p_value = test$p.value
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  select(Genus, Variable, rho, p_value, p_adj)

write.csv(cor_results, file.path(out_tab, "TABLE_Spearman_Genus_vs_Env.csv"), row.names = FALSE)

sig_cor <- cor_results %>% filter(p_adj <= 0.05, abs(rho) >= 0.4)
write.csv(sig_cor, file.path(out_tab, "TABLE_Significant_Genus_vs_Env_FDR0.05_rho0.4.csv"), row.names = FALSE)

df_plot <- sig_cor
if(nrow(df_plot) == 0){
  message("No correlations pass FDR<=0.05 & |rho|>=0.4; plotting top |rho| genera per variable instead.")
  df_plot <- cor_results %>% group_by(Variable) %>% slice_max(order_by = abs(rho), n = 20) %>% ungroup()
}

df_plot <- df_plot %>%
  mutate(
    Genus = fct_reorder(Genus, abs(rho), .fun = max),
    Variable = factor(Variable, levels = unique(Variable))
  )

p_cor <- ggplot(df_plot, aes(x = Variable, y = Genus, fill = rho)) +
  geom_tile(color = "grey80", linewidth = 0.3) +
  scale_fill_gradient2(midpoint = 0, labels = scales::number_format(accuracy = 0.01)) +
  labs(title = "Genus ↔ Metals / Organic Matter (Spearman)",
       subtitle = "rho shown; full table exported with FDR-adjusted p-values",
       x = "Environmental variable", y = "Genus", fill = "rho") +
  theme(panel.grid = element_blank())

ggsave(file.path(out_fig, "FIG_16S_Genus_Env_CorrelationHeatmap.png"),
       p_cor, width = 9, height = 10, dpi = 300)
ggsave(file.path(out_fig, "FIG_16S_Genus_Env_CorrelationHeatmap.pdf"),
       p_cor, width = 9, height = 10)

### =========================
### 13) REPLICATE SITE SUMMARY TABLE (ENV MEANS + GENUS MEAN REL)
### =========================
site_env <- meta_df %>%
  group_by(site_abbrev) %>%
  summarize(
    zone_mode = names(sort(table(as.character(zone)), decreasing = TRUE))[1],
    across(all_of(env_cols), ~mean(as.numeric(.x), na.rm = TRUE), .names = "mean_{.col}"),
    .groups = "drop"
  )

site_genus_wide <- site_genus %>%
  pivot_wider(names_from = Genus, values_from = mean_rel, values_fill = 0)

site_summary <- site_env %>%
  left_join(site_genus_wide, by = "site_abbrev")

write.csv(site_summary, file.path(out_tab, "TABLE_SiteAbbrevSummary_EnvPlusGenusMeanRel.csv"), row.names = FALSE)

### =========================
### DONE
### =========================
message("✅ Pipeline complete (replicates separate).")
message("Figures saved to: ", out_fig)
message("Tables  saved to: ", out_tab)
message("Sample mapping table saved to: ", file.path(out_tab, "TABLE_SampleName_Mapping_ASV_to_Metadata.csv"))
###############################################
