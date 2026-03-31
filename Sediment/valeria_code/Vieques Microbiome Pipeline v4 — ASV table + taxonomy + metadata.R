# ===============================================================
# Vieques Microbiome Pipeline v4 — ASV table + taxonomy + metadata
# Author: auto-generated for Valeria
# Inputs (edit the paths if needed):
#   - "asv_table.csv"
#   - "asv_id_map_with_taxonomy.csv"                (ID column = 'asv')
#   - "Meta data  - Copy of sediment toxicity .csv" (contains SampleID)
#   - "CHEM_SED.xlsx"
# Outputs: exports/ folder with CSV tables, plots (PNG), and diagnostics.
# ===============================================================

# ---------------------------
# 0) Packages & helpers
# ---------------------------
need <- c(
  "tidyverse","data.table","janitor","readxl","stringr","glue",
  "vegan","phyloseq","ape","matrixStats","broom",
  "ggpubr","FSA","rstatix","pheatmap","scales","magrittr","tibble"
)
for(p in need){
  if(!requireNamespace(p, quietly = TRUE)){
    install.packages(p, repos="https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}
set.seed(42)
theme_set(theme_bw(base_size = 12))
dir.create("exports", showWarnings = FALSE)

# Helpers
stdz <- function(x) as.numeric(scale(x))
lz   <- function(x) stdz(log10(x + 1))
first_existing_col <- function(df, candidates) {
  cand <- intersect(candidates, names(df))
  if (length(cand)) cand[1] else NULL
}
coalesce_across <- function(df, candidates) {
  cand <- intersect(candidates, names(df))
  if (!length(cand)) return(rep(NA_character_, nrow(df)))
  out <- df[[cand[1]]]
  if (length(cand) > 1) {
    for (i in cand[-1]) {
      out <- dplyr::coalesce(dplyr::na_if(out, ""), dplyr::na_if(df[[i]], ""))
    }
  }
  out
}
clean_id <- function(x){
  x %>% as.character() %>% trimws() %>% str_replace_all("\\s+", "_")
}

# ---------------------------
# 1) User paths (EDIT if needed)
# ---------------------------
asv_path  <- "asv_table.csv"
tax_path  <- "asv_id_map_with_taxonomy.csv"
meta_path <- "Meta data  - Copy of sediment toxicity .csv"  # spaces intentional
chem_path <- "CHEM_SED.xlsx"

stopifnot(file.exists(asv_path), file.exists(tax_path),
          file.exists(meta_path), file.exists(chem_path))
message("✅ Inputs found.")

# ---------------------------
# 2) Chemistry → indices (PC1, z-mean, tiers)
# ---------------------------
chem_raw <- readxl::read_excel(chem_path) |> janitor::clean_names()

chem <- chem_raw |>
  rename(
    site = any_of(c("site","site_name")),
    zone = any_of(c("zone")),
    site_abbrev = any_of(c("site_abbreviation","site_abbrev","site_abbr","site__abbreviation",
                           "site__abbr","site_abbreviation_")),
    cu = any_of(c("cu_mg_kg","cu_mgkg","cu","cu_mg_k_g")),
    ni = any_of(c("ni_mg_kg","ni_mgkg","ni")),
    zn = any_of(c("zn_mg_kg","zn_mgkg","zn")),
    co = any_of(c("co_mg_kg","co_mgkg","co"))
  )

# fallback metal grabs if weird names
nm <- names(chem_raw)
grab_metal <- function(keys){
  cand <- nm[str_detect(nm, keys)]
  if(length(cand)) chem_raw[[cand[1]]] else NA_real_
}
if(!"cu" %in% names(chem)) chem$cu <- grab_metal("cu.*mg.*kg|\\bcu\\b")
if(!"ni" %in% names(chem)) chem$ni <- grab_metal("ni.*mg.*kg|\\bni\\b")
if(!"zn" %in% names(chem)) chem$zn <- grab_metal("zn.*mg.*kg|\\bzn\\b")
if(!"co" %in% names(chem)) chem$co <- grab_metal("co.*mg.*kg|\\bco\\b")

# site_abbrev derive if needed
if(is.null(chem$site_abbrev)){
  chem$site_abbrev <- coalesce_across(chem_raw,
                                      c("site_abbreviation","site_abbrev","abbrev","abbr","site_abbr","site__abbreviation"))
  if(is.null(chem$site_abbrev) && "site" %in% names(chem)) {
    chem$site_abbrev <- str_replace_all(toupper(substr(chem$site,1,4)), "[^A-Z0-9]", "")
  }
}

chem <- chem |>
  mutate(across(c(cu,ni,zn,co), as.numeric)) |>
  select(site, zone, site_abbrev, cu, ni, zn, co) |>
  distinct()

metals_trans <- chem |>
  select(cu,ni,zn,co) |>
  mutate(across(everything(), lz))

chem$contam_index_zmean <- rowMeans(metals_trans, na.rm = TRUE)
pc <- prcomp(metals_trans |> mutate(across(everything(), ~replace_na(.,0))),
             center = TRUE, scale. = FALSE)
chem$contam_pc1 <- pc$x[,1]

q <- quantile(chem$contam_pc1, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE) |> unique()
if (length(q) < 4) { rng <- range(chem$contam_pc1, na.rm = TRUE); q <- seq(rng[1], rng[2], length.out = 4) }
chem$contam_tier <- cut(chem$contam_pc1, breaks = q, include.lowest = TRUE,
                        labels = c("Low","Medium","High"))

readr::write_csv(chem, "exports/chem_clean_with_index.csv")

# ---------------------------
# 3) Metadata (load first; normalize IDs)
# ---------------------------
meta <- suppressMessages(readr::read_csv(meta_path, guess_max = 100000)) %>% janitor::clean_names()
sample_col <- first_existing_col(meta, c("sampleid","sample_id","sample","id","sample_name"))
if (is.null(sample_col)) stop("SampleID column not found in metadata. Available: ", paste(names(meta), collapse=", "))
names(meta)[names(meta)==sample_col] <- "SampleID"
meta$SampleID <- clean_id(meta$SampleID)

site_candidates <- c("site_abbrev","site_abbreviation","site_abbr","abbrev","abbr","site","site_name")
zone_col <- first_existing_col(meta, c("zone","Zone"))

meta <- meta %>%
  mutate(
    site_abbrev = coalesce_across(., site_candidates),
    site_abbrev = toupper(trimws(as.character(site_abbrev))),
    zone = if (!is.null(zone_col)) .[[zone_col]] else NA_character_
  )

# ---------------------------
# 4) ASV table + Taxonomy → phyloseq
# ---------------------------
# ASV table: auto-detect orientation and transpose if samples are rows
asv_raw <- suppressMessages(readr::read_csv(asv_path))
first_col <- names(asv_raw)[1]
looks_like_taxa   <- grepl("^(ASV|Feature|OTU|Seq|Sequence)", first_col, ignore.case = TRUE)
looks_like_sample <- grepl("^(Sample|SampleID|ID|Run|Library)", first_col, ignore.case = TRUE)

if (looks_like_sample && !looks_like_taxa) {
  message("Detected samples in rows; transposing ASV table so taxa are rows.")
  mat <- as.matrix(asv_raw[,-1])
  rownames(mat) <- clean_id(asv_raw[[1]])
  mat[is.na(mat)] <- 0
  mat <- t(mat)
  asv_df <- as.data.frame(mat) %>% tibble::rownames_to_column("ASV_ID")
} else {
  asv_df <- asv_raw
  if (!grepl("ASV|Feature|OTU|Seq", first_col, ignore.case = TRUE)) names(asv_df)[1] <- "ASV_ID"
  asv_df[[1]] <- clean_id(asv_df[[1]])
}
# numeric counts
for (j in 2:ncol(asv_df)) asv_df[[j]] <- suppressWarnings(as.numeric(asv_df[[j]]))
asv_df[is.na(asv_df)] <- 0

# Taxonomy: your file key column is 'asv'
tax_raw <- suppressMessages(readr::read_csv(tax_path)) %>% janitor::clean_names()
if (!"asv" %in% names(tax_raw)) stop("Taxonomy file must contain column 'asv'. Found: ", paste(names(tax_raw), collapse=", "))
tax_raw$asv <- clean_id(tax_raw$asv)

tax_df <- tax_raw %>%
  rename(ASV_ID = asv) %>%
  select(ASV_ID, any_of(c("kingdom","phylum","class","order","family","genus","species"))) %>%
  distinct()

# Build matrices
asv_ids <- asv_df$ASV_ID
otu_mat <- as.matrix(asv_df[,-1]); rownames(otu_mat) <- asv_ids

tax_mat <- tax_df %>%
  right_join(tibble(ASV_ID = asv_ids), by = "ASV_ID") %>%
  select(-ASV_ID) %>% as.matrix()
rownames(tax_mat) <- asv_ids

# Harmonize samples with metadata
colnames(otu_mat) <- clean_id(colnames(otu_mat))
meta$SampleID     <- clean_id(meta$SampleID)
common_samples <- intersect(colnames(otu_mat), meta$SampleID)
message(glue::glue("Samples in ASV table:   {ncol(otu_mat)}"))
message(glue::glue("Samples in metadata:    {nrow(meta)}"))
message(glue::glue("Overlap (by SampleID):  {length(common_samples)}"))
if (length(common_samples) < 2) stop("Fewer than 2 overlapping samples between ASV table and metadata after ID cleaning.")

# Subset & order
otu_mat <- otu_mat[, common_samples, drop=FALSE]
meta2   <- meta %>% filter(SampleID %in% common_samples) %>% distinct(SampleID, .keep_all = TRUE) %>% as.data.frame()
rownames(meta2) <- meta2$SampleID
meta2$zone <- as.factor(meta2$zone)

# phyloseq object
OTU  <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX  <- tax_table(tax_mat)
META <- sample_data(meta2)
ps   <- phyloseq(OTU, TAX, META)

# prune zeros
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)
stopifnot(nsamples(ps) > 1, ntaxa(ps) > 1)

# ---------------------------
# 5) Join chemistry onto samples (by site_abbrev)
# ---------------------------
# ---------------------------
# 5) Join chemistry onto samples (by site_abbrev)
# ---------------------------

# 5a) Start from a plain data.frame, not sample_data()
sd_df <- as(sample_data(ps), "data.frame")
sd_df$SampleID <- rownames(sd_df)

# Standardize key in chemistry
chem_merge <- chem %>%
  mutate(site_abbrev = toupper(trimws(as.character(site_abbrev))))

# 5b) Left-join chem → samples by site_abbrev
sd2_df <- sd_df %>%
  dplyr::left_join(chem_merge, by = "site_abbrev")

# Report join success
message(glue::glue(
  "Joined chemistry for {sum(!is.na(sd2_df$contam_pc1))}/{nrow(sd2_df)} samples (by site_abbrev)."
))

# ---- 5c’) Robustly synthesize a single `zone` column, then write CSV and assign ----

# Start with a blank character vector of the right length
sd2_df$zone <- rep(NA_character_, nrow(sd2_df))

# If zone.x exists (from metadata), prefer it when non-empty
if ("zone.x" %in% names(sd2_df)) {
  z1 <- as.character(sd2_df$zone.x)
  keep1 <- !is.na(z1) & nzchar(z1)
  sd2_df$zone[keep1] <- z1[keep1]
}

# If zone.y exists (from chemistry), fill any remaining NAs/empties
if ("zone.y" %in% names(sd2_df)) {
  z2 <- as.character(sd2_df$zone.y)   # coerce to character to avoid type clashes
  need2 <- is.na(sd2_df$zone) | !nzchar(sd2_df$zone)
  sd2_df$zone[need2] <- z2[need2]
}

# Make zone a factor for downstream models
sd2_df$zone <- factor(sd2_df$zone)

# Tidy up: drop the suffix columns if present
sd2_df <- sd2_df %>% dplyr::select(-dplyr::any_of(c("zone.x","zone.y")))

# Diagnostics CSV
readr::write_csv(
  sd2_df %>% dplyr::transmute(SampleID, site_abbrev, zone, has_chem = !is.na(contam_pc1)),
  "exports/join_diagnostics_samples_vs_chem.csv"
)

# Assign back to phyloseq
sd2_ps <- sd2_df
rownames(sd2_ps) <- sd2_ps$SampleID
sample_data(ps) <- phyloseq::sample_data(sd2_ps)

# ---- 5c’) Robustly synthesize a single `zone` column, then write CSV and assign ----

# Start with a blank character vector of the right length
sd2_df$zone <- rep(NA_character_, nrow(sd2_df))

# If zone.x exists (from metadata), prefer it when non-empty
if ("zone.x" %in% names(sd2_df)) {
  z1 <- as.character(sd2_df$zone.x)
  keep1 <- !is.na(z1) & nzchar(z1)
  sd2_df$zone[keep1] <- z1[keep1]
}

# If zone.y exists (from chemistry), fill any remaining NAs/empties
if ("zone.y" %in% names(sd2_df)) {
  z2 <- as.character(sd2_df$zone.y)   # coerce to character to avoid type clashes
  need2 <- is.na(sd2_df$zone) | !nzchar(sd2_df$zone)
  sd2_df$zone[need2] <- z2[need2]
}

# Make zone a factor for downstream models
sd2_df$zone <- factor(sd2_df$zone)

# Tidy up: drop the suffix columns if present
sd2_df <- sd2_df %>% dplyr::select(-dplyr::any_of(c("zone.x","zone.y")))

# Diagnostics CSV
readr::write_csv(
  sd2_df %>% dplyr::transmute(SampleID, site_abbrev, zone, has_chem = !is.na(contam_pc1)),
  "exports/join_diagnostics_samples_vs_chem.csv"
)

# Assign back to phyloseq
sd2_ps <- sd2_df
rownames(sd2_ps) <- sd2_ps$SampleID
sample_data(ps) <- phyloseq::sample_data(sd2_ps)

message("zone source → ",
        if ("zone.x" %in% names(sd2_df)) "metadata " else "",
        if ("zone.y" %in% names(sd2_df)) "+ chemistry" else "")









# ---------------------------
# 6) Alpha diversity
# ---------------------------
# ---------------------------
# 6) Alpha diversity
# ---------------------------
# ---------------------------
# 6) Alpha diversity  (robust SampleID handling)
# ---------------------------

# Get plain data.frame of sample metadata
sd_alpha_df <- as(phyloseq::sample_data(ps), "data.frame")

# Add a fallback rowname column (won't collide with existing SampleID)
sd_alpha_df$SampleID_rowname <- rownames(sd_alpha_df)

# Use existing SampleID if present; otherwise use rownames
if ("SampleID" %in% names(sd_alpha_df)) {
  sd_alpha_df$SampleID <- clean_id(sd_alpha_df$SampleID)
} else {
  sd_alpha_df$SampleID <- clean_id(sd_alpha_df$SampleID_rowname)
}

# Keep only the fields we need (and drop the helper column)
sd_alpha_df <- sd_alpha_df |>
  dplyr::select(
    SampleID,
    dplyr::any_of(c("site_abbrev","zone","contam_pc1","contam_index_zmean","contam_tier"))
  )

# Compute alpha and join metadata
alpha_df <- phyloseq::estimate_richness(ps, measures = c("Observed","Chao1","Shannon")) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::left_join(sd_alpha_df, by = "SampleID")

readr::write_csv(alpha_df, "exports/alpha_metrics_merged.csv")

# Rows complete for each analysis
alpha_cc_tier <- alpha_df |> dplyr::filter(!is.na(contam_tier))
alpha_cc_pc1  <- alpha_df |> dplyr::filter(!is.na(contam_pc1))

# Plots
g_alpha_box <- function(metric){
  ggplot(alpha_cc_tier, aes(x = contam_tier, y = .data[[metric]], fill = contam_tier)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_jitter(width = 0.08, height = 0, alpha = 0.7, size = 1) +
    labs(x = "Contamination tier (PC1 tertiles)", y = metric,
         title = glue::glue("{metric} by contamination tier")) +
    guides(fill="none")
}
ggsave("exports/alpha_Shannon_by_tier.png", g_alpha_box("Shannon"), width=6, height=4, dpi=300)
ggsave("exports/alpha_Chao1_by_tier.png",   g_alpha_box("Chao1"),   width=6, height=4, dpi=300)
ggsave("exports/alpha_Observed_by_tier.png",g_alpha_box("Observed"),width=6, height=4, dpi=300)

# Stats
alpha_tests <- dplyr::bind_rows(
  rstatix::kruskal_test(alpha_cc_tier, Shannon ~ contam_tier) |> dplyr::mutate(metric="Shannon"),
  rstatix::kruskal_test(alpha_cc_tier, Chao1   ~ contam_tier) |> dplyr::mutate(metric="Chao1"),
  rstatix::kruskal_test(alpha_cc_tier, Observed~ contam_tier) |> dplyr::mutate(metric="Observed")
)
readr::write_csv(alpha_tests, "exports/alpha_kruskal_by_tier.csv")

dunn_tbl <- list(
  Shannon = FSA::dunnTest(Shannon ~ contam_tier, data = alpha_cc_tier, method = "bh")$res,
  Chao1   = FSA::dunnTest(Chao1   ~ contam_tier, data = alpha_cc_tier, method = "bh")$res,
  Observed= FSA::dunnTest(Observed~ contam_tier, data = alpha_cc_tier, method = "bh")$res
) |> purrr::imap(~dplyr::mutate(.x, metric = .y)) |> dplyr::bind_rows()
readr::write_csv(dunn_tbl, "exports/alpha_dunn_by_tier.csv")

# Linear models vs PC1
# ---------------------------
# Linear models vs PC1 (robust to single-level zone)
# ---------------------------

safe_lm_alpha <- function(metric){
  d <- alpha_cc_pc1 |>
    dplyr::select(y = dplyr::all_of(metric), contam_pc1, zone) |>
    tidyr::drop_na(y, contam_pc1)
  
  # ensure factor, drop empty levels
  d$zone <- droplevels(factor(d$zone))
  has_zone <- nlevels(d$zone) >= 2
  has_var  <- stats::var(d$contam_pc1, na.rm = TRUE) > 0
  enough_n <- nrow(d) >= 3
  
  if (!enough_n || !has_var) {
    return(
      tibble::tibble(
        term      = c("(Intercept)", "contam_pc1"),
        estimate  = NA_real_, std.error = NA_real_,
        statistic = NA_real_, p.value  = NA_real_,
        metric    = metric,
        model     = if (has_zone) "y ~ contam_pc1 + zone" else "y ~ contam_pc1",
        note      = dplyr::case_when(
          !enough_n ~ "Skipped: <3 complete rows after NA drop",
          !has_var  ~ "Skipped: contam_pc1 has zero variance",
          TRUE ~ "Skipped"
        )
      )
    )
  }
  
  form <- if (has_zone) y ~ contam_pc1 + zone else y ~ contam_pc1
  fit  <- stats::lm(form, data = d)
  
  broom::tidy(fit) |>
    dplyr::mutate(
      metric = metric,
      model  = if (has_zone) "y ~ contam_pc1 + zone" else "y ~ contam_pc1",
      note   = dplyr::if_else(has_zone, NA_character_, "zone had <2 levels; dropped")
    )
}

lm_out <- dplyr::bind_rows(
  safe_lm_alpha("Shannon"),
  safe_lm_alpha("Chao1"),
  safe_lm_alpha("Observed")
)

readr::write_csv(lm_out, "exports/alpha_lm_vs_pc1_zone.csv")

# Scatter plots
g_alpha_scatter <- function(metric){
  ggplot(alpha_cc_pc1, aes(x = contam_pc1, y = .data[[metric]], color = zone)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(x = "Contamination gradient (PC1)", y = metric,
         title = glue::glue("{metric} ~ contamination (PC1)")) +
    theme(legend.position = "bottom")
}
ggsave("exports/alpha_Shannon_vs_pc1.png", g_alpha_scatter("Shannon"), width=6, height=4, dpi=300)
ggsave("exports/alpha_Chao1_vs_pc1.png",   g_alpha_scatter("Chao1"),   width=6, height=4, dpi=300)


# ---------------------------
# 7) Beta diversity: PERMANOVA, betadisper, dbRDA, envfit
# ---------------------------
# ---------------------------
# 7) Beta diversity: PERMANOVA, betadisper, dbRDA, envfit (robust)
# ---------------------------

# ---- Helpers (scoped to this section) ----
run_adonis_safe <- function(dist_mat, data_df, formula_str, outfile){
  ok_n <- nrow(data_df) >= 3 && nrow(dist_mat) >= 3
  if (!ok_n) {
    writeLines(glue::glue("PERMANOVA skipped: need >=3 samples (have {nrow(data_df)})."),
               outfile)
    return(invisible(NULL))
  }
  fit <- tryCatch(
    vegan::adonis2(as.formula(formula_str), data = data_df, permutations = 9999, by = "margin"),
    error = function(e) e
  )
  sink(outfile)
  cat("PERMANOVA (Bray)\n")
  cat("Model: ", formula_str, "\n", sep = "")
  if (inherits(fit, "error")) {
    cat("Skipped with error: ", conditionMessage(fit), "\n", sep = "")
  } else {
    print(fit)
  }
  sink()
}

run_betadisper_safe <- function(dist_mat, groups, outfile, label="group"){
  df <- data.frame(g = droplevels(groups))
  tbl <- table(df$g)
  enough <- length(tbl[tbl >= 2]) >= 2 && nrow(df) >= 4 && nrow(dist_mat) >= 4
  if (!enough) {
    msg <- paste0(
      "Betadisper skipped (", label, "). Need >=2 groups with >=2 samples.\n",
      "Group sizes: ", paste(paste(names(tbl), tbl, sep=":"), collapse=", "), "\n"
    )
    writeLines(msg, outfile)
    return(invisible(NULL))
  }
  bd <- tryCatch(betadisper(as.dist(dist_mat), df$g), error = function(e) e)
  sink(outfile)
  cat("Betadisper by ", label, "\n", sep = "")
  if (inherits(bd, "error")) {
    cat("Skipped with error: ", conditionMessage(bd), "\n", sep = "")
  } else {
    print(anova(bd))
  }
  sink()
}

# ---- Distances & ordination ----
ps_rel    <- transform_sample_counts(ps, function(x) if (sum(x)==0) x else x / sum(x))
dist_bray <- phyloseq::distance(ps_rel, method = "bray")
ord_pcoa  <- ordinate(ps_rel, method = "PCoA", distance = dist_bray)

# ---- Metadata frame & types ----
mdf <- as(sample_data(ps_rel), "data.frame")
mdf$SampleID     <- rownames(mdf)
mdf$zone         <- droplevels(factor(mdf$zone))
mdf$site_abbrev  <- droplevels(factor(mdf$site_abbrev))
mdf$contam_tier  <- droplevels(factor(mdf$contam_tier))
mdf$contam_pc1   <- suppressWarnings(as.numeric(mdf$contam_pc1))

# ---------- Overall PERMANOVA (chooses best feasible model) ----------
keep_all  <- !is.na(mdf$contam_pc1)
mdf_all   <- mdf[keep_all, , drop = FALSE]
dist_all  <- as.matrix(dist_bray)[keep_all, keep_all, drop = FALSE]

has_zone  <- nlevels(mdf_all$zone)        >= 2
has_site  <- nlevels(mdf_all$site_abbrev) >= 2
has_var   <- is.finite(var(mdf_all$contam_pc1, na.rm=TRUE)) && var(mdf_all$contam_pc1, na.rm=TRUE) > 0

form_str <- if (has_zone && has_site) {
  "as.dist(dist_all) ~ contam_pc1 + zone + site_abbrev"
} else if (has_zone) {
  "as.dist(dist_all) ~ contam_pc1 + zone"
} else if (has_site) {
  "as.dist(dist_all) ~ contam_pc1 + site_abbrev"
} else {
  "as.dist(dist_all) ~ contam_pc1"
}

if (!has_var) {
  writeLines("PERMANOVA skipped: contam_pc1 variance == 0.",
             "exports/permanova_overall.txt")
} else {
  run_adonis_safe(dist_all, mdf_all, form_str, "exports/permanova_overall.txt")
}

# ---------- Dispersion tests ----------
run_betadisper_safe(dist_all, mdf_all$zone,        "exports/betadisper_by_zone.txt",     label="zone")
run_betadisper_safe(dist_all, mdf_all$site_abbrev, "exports/betadisper_by_location.txt", label="site_abbrev")

# ---------- dbRDA (capscale) ----------
if (nrow(mdf_all) >= 3 && nrow(dist_all) >= 3 && has_var) {
  form_cap <- if (has_zone && has_site) {
    "as.dist(dist_all) ~ contam_pc1 + zone + site_abbrev"
  } else if (has_zone) {
    "as.dist(dist_all) ~ contam_pc1 + zone"
  } else if (has_site) {
    "as.dist(dist_all) ~ contam_pc1 + site_abbrev"
  } else {
    "as.dist(dist_all) ~ contam_pc1"
  }
  cap <- tryCatch(capscale(as.formula(form_cap), data = mdf_all, add = TRUE),
                  error = function(e) e)
  if (!inherits(cap, "error")) {
    saveRDS(cap, "exports/capscale_bray_model.rds")
    png("exports/ordination_dbRDA_bray.png", width=1600, height=1200, res=200)
    plot(cap, display = c("sites","bp"),
         main = paste0("dbRDA (Bray) ~ ", gsub("as.dist\\(dist_all\\) ~ ", "", form_cap)))
    dev.off()
  } else {
    writeLines(paste0("dbRDA skipped: ", conditionMessage(cap)),
               "exports/ordination_dbRDA_bray.txt")
  }
} else {
  writeLines("dbRDA skipped: insufficient usable samples or zero variance in contam_pc1.",
             "exports/ordination_dbRDA_bray.txt")
}

# ---------- envfit on unconstrained PCoA ----------
env_vars <- mdf %>%
  dplyr::select(cu, ni, zn, co, contam_pc1) %>%
  dplyr::mutate(dplyr::across(everything(), ~ suppressWarnings(as.numeric(.x))))
ef <- tryCatch(vegan::envfit(ord_pcoa, env_vars, permutations = 9999),
               error = function(e) e)
sink("exports/envfit_pcoa_bray.txt")
if (inherits(ef, "error")) {
  cat("envfit skipped: ", conditionMessage(ef), "\n", sep = "")
} else {
  print(ef)
}
sink()

# ---------- Per-location PERMANOVA (subset within each site_abbrev) ----------
dir.create("exports/per_location", showWarnings = FALSE)
sites <- levels(mdf$site_abbrev)
per_site_summary <- list()

for (s in sites) {
  idx <- which(mdf$site_abbrev == s & !is.na(mdf$contam_pc1))
  if (length(idx) < 3) {
    writeLines(glue::glue("Site {s}: skipped (n={length(idx)})."),
               file.path("exports/per_location", paste0("permanova_site_", s, ".txt")))
    next
  }
  mdf_s  <- mdf[idx, , drop = FALSE]
  dist_s <- as.matrix(dist_bray)[idx, idx, drop = FALSE]
  
  has_zone_s <- nlevels(droplevels(mdf_s$zone)) >= 2
  has_var_s  <- is.finite(var(mdf_s$contam_pc1, na.rm=TRUE)) && var(mdf_s$contam_pc1, na.rm=TRUE) > 0
  
  if (!has_var_s) {
    writeLines(glue::glue("Site {s}: PERMANOVA skipped (contam_pc1 variance == 0)."),
               file.path("exports/per_location", paste0("permanova_site_", s, ".txt")))
    next
  }
  
  form_s   <- if (has_zone_s) "as.dist(dist_s) ~ contam_pc1 + zone" else "as.dist(dist_s) ~ contam_pc1"
  out_path <- file.path("exports/per_location", paste0("permanova_site_", s, ".txt"))
  run_adonis_safe(dist_s, mdf_s, form_s, out_path)
  
  per_site_summary[[s]] <- data.frame(
    site_abbrev = s,
    n = nrow(mdf_s),
    zone_levels = nlevels(droplevels(mdf_s$zone)),
    used_model  = gsub("as.dist\\(dist_s\\) ~ ", "", form_s),
    stringsAsFactors = FALSE
  )
}
if (length(per_site_summary)) {
  per_site_summary <- dplyr::bind_rows(per_site_summary)
  readr::write_csv(per_site_summary, "exports/per_location/per_site_summary.csv")
}

# ---------- PCoA plots (by zone / by location / faceted) ----------
scores <- as.data.frame(ord_pcoa$vectors)
scores$SampleID <- rownames(scores)
scores <- scores %>%
  dplyr::left_join(mdf %>% dplyr::select(SampleID, zone, site_abbrev, contam_tier), by="SampleID")

pc1 <- round(ord_pcoa$values$Relative_eig[1]*100,1)
pc2 <- round(ord_pcoa$values$Relative_eig[2]*100,1)

# Color by zone
p_zone <- ggplot(scores, aes(Axis.1, Axis.2, color = zone)) +
  geom_point(size = 2, alpha = 0.9) +
  stat_ellipse(level = 0.68, linetype = 2, linewidth = 0.7, alpha = 0.6, na.rm = TRUE) +
  labs(x = glue::glue("PCoA1 ({pc1}%)"), y = glue::glue("PCoA2 ({pc2}%)"),
       title = "PCoA (Bray-Curtis) colored by zone") +
  theme(legend.position = "bottom")
ggsave("exports/beta_pcoa_by_zone.png", p_zone, width=7, height=5, dpi=300)

# Color by location
p_loc <- ggplot(scores, aes(Axis.1, Axis.2, color = site_abbrev)) +
  geom_point(size = 2, alpha = 0.9) +
  stat_ellipse(level = 0.68, linetype = 2, linewidth = 0.7, alpha = 0.6, na.rm = TRUE) +
  labs(x = glue::glue("PCoA1 ({pc1}%)"), y = glue::glue("PCoA2 ({pc2}%)"),
       title = "PCoA (Bray-Curtis) colored by location") +
  theme(legend.position = "bottom")
ggsave("exports/beta_pcoa_by_location.png", p_loc, width=7, height=5, dpi=300)

# Facet by zone
p_facet_zone <- ggplot(scores, aes(Axis.1, Axis.2, color = site_abbrev)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(x = glue::glue("PCoA1 ({pc1}%)"), y = glue::glue("PCoA2 ({pc2}%)"),
       title = "PCoA (Bray-Curtis) faceted by zone") +
  theme(legend.position = "bottom") +
  facet_wrap(~ zone, drop = TRUE)
ggsave("exports/beta_pcoa_facet_by_zone.png", p_facet_zone, width=9, height=6, dpi=300)

# Facet by location
p_facet_loc <- ggplot(scores, aes(Axis.1, Axis.2, color = zone)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(x = glue::glue("PCoA1 ({pc1}%)"), y = glue::glue("PCoA2 ({pc2}%)"),
       title = "PCoA (Bray-Curtis) faceted by location") +
  theme(legend.position = "bottom") +
  facet_wrap(~ site_abbrev, drop = TRUE)
ggsave("exports/beta_pcoa_facet_by_location.png", p_facet_loc, width=10, height=7, dpi=300)

# Export scores
readr::write_csv(scores, "exports/pcoa_scores_with_zone_and_location.csv")


# ---------------------------
# 8) Taxa ~ contaminants (correlations + optional ANCOM-BC2)
# ---------------------------
# ---------------------------
# 8) Taxa ~ contaminants (correlations + optional ANCOM-BC2) — ROBUST
# ---------------------------

# Helper: pick deepest rank with enough non-NA/non-empty labels
pick_rank <- function(ps, candidates = c("Genus","genus","Family","family","Order","order","Class","class","Phylum","phylum","Kingdom","kingdom")) {
  if (is.null(tax_table(ps, errorIfNULL = FALSE))) return(NA_character_)
  tt <- as.data.frame(tax_table(ps))
  for (rk in candidates) {
    if (rk %in% colnames(tt)) {
      vals <- as.character(tt[[rk]])
      nn <- sum(!is.na(vals) & nzchar(vals) & vals != "NA")
      if (nn >= 2) return(rk)
    }
  }
  NA_character_
}

collapse_safe <- function(ps, rk) {
  if (is.na(rk)) return(list(obj = ps, used = "ASV (no collapse; no usable rank)"))
  # keep NA-labeled taxa together instead of dropping them
  tmp <- try(suppressWarnings(phyloseq::tax_glom(ps, taxrank = rk, NArm = FALSE)), silent = TRUE)
  if (inherits(tmp, "try-error") || ntaxa(tmp) == 0) {
    list(obj = ps, used = paste0("ASV (fallback; collapse at ", rk, " failed)"))
  } else {
    list(obj = tmp, used = rk)
  }
}

rk <- pick_rank(ps)
col_res <- collapse_safe(ps, rk)
ps_coll <- col_res$obj
used_rank <- col_res$used
message("Taxa collapse level used: ", used_rank)

# Relative abundance for correlations
ps_coll_rel <- transform_sample_counts(ps_coll, function(x) if (sum(x)==0) x else x / sum(x))

# Pull tables
otu_df <- as.data.frame(otu_table(ps_coll_rel))
tax_df <- if (!is.null(tax_table(ps_coll_rel, errorIfNULL = FALSE))) {
  as.data.frame(tax_table(ps_coll_rel)) %>% tibble::rownames_to_column("Feature")
} else {
  tibble::tibble(Feature = rownames(otu_df))
}

# Long table: SampleID x Feature x rel_abund
long_df <- otu_df %>%
  tibble::rownames_to_column("Feature") %>%
  tidyr::pivot_longer(-Feature, names_to = "SampleID", values_to = "rel_abund")

# Choose a display name column for taxa
tax_disp_col <- intersect(c("Genus","genus","Family","family","Order","order","Class","class","Phylum","phylum","Kingdom","kingdom"), colnames(tax_df))
if (length(tax_disp_col) == 0) tax_disp_col <- "Feature"
disp_name <- if (length(tax_disp_col) > 0) tax_disp_col[1] else "Feature"

long_df <- long_df %>%
  dplyr::left_join(tax_df %>% dplyr::select(Feature, !!disp_name), by = "Feature") %>%
  dplyr::rename(Taxon = !!disp_name)

# Numerical env variables from the Section 7 metadata (mdf)
stopifnot(exists("mdf"))
num_vars <- c("cu","ni","zn","co","contam_pc1")
mdf_num <- mdf %>%
  dplyr::transmute(SampleID = rownames(.),
                   cu = suppressWarnings(as.numeric(cu)),
                   ni = suppressWarnings(as.numeric(ni)),
                   zn = suppressWarnings(as.numeric(zn)),
                   co = suppressWarnings(as.numeric(co)),
                   contam_pc1 = suppressWarnings(as.numeric(contam_pc1)))

# Spearman correlations Taxon vs metals/PC1 (per taxon), BH FDR
safe_cor_test <- function(x, y) {
  # Handle constant vectors or all-NA gracefully
  if (all(is.na(x)) || all(is.na(y))) return(c(rho = NA_real_, p = NA_real_))
  if (stats::var(x, na.rm = TRUE) == 0 || stats::var(y, na.rm = TRUE) == 0) return(c(rho = NA_real_, p = NA_real_))
  ct <- try(stats::cor.test(x, y, method = "spearman"), silent = TRUE)
  if (inherits(ct, "try-error")) return(c(rho = NA_real_, p = NA_real_))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

gen_cor <- long_df %>%
  dplyr::filter(!is.na(Taxon) & nzchar(Taxon)) %>%
  dplyr::left_join(mdf_num, by = "SampleID") %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarize(
    # rho
    rho_cu = safe_cor_test(rel_abund, cu)["rho"],
    rho_ni = safe_cor_test(rel_abund, ni)["rho"],
    rho_zn = safe_cor_test(rel_abund, zn)["rho"],
    rho_co = safe_cor_test(rel_abund, co)["rho"],
    rho_contam_pc1 = safe_cor_test(rel_abund, contam_pc1)["rho"],
    # p
    p_cu = safe_cor_test(rel_abund, cu)["p"],
    p_ni = safe_cor_test(rel_abund, ni)["p"],
    p_zn = safe_cor_test(rel_abund, zn)["p"],
    p_co = safe_cor_test(rel_abund, co)["p"],
    p_contam_pc1 = safe_cor_test(rel_abund, contam_pc1)["p"],
    .groups = "drop"
  ) %>%
  # FDR for each p-column
  dplyr::mutate(
    q_cu = p.adjust(p_cu, method = "BH"),
    q_ni = p.adjust(p_ni, method = "BH"),
    q_zn = p.adjust(p_zn, method = "BH"),
    q_co = p.adjust(p_co, method = "BH"),
    q_contam_pc1 = p.adjust(p_contam_pc1, method = "BH"),
    collapse_level = used_rank
  ) %>%
  # For safety, coerce to numeric
  dplyr::mutate(dplyr::across(dplyr::starts_with(c("rho_","p_","q_")), ~ as.numeric(.)))

readr::write_csv(gen_cor, "exports/genus_spearman_correlations_vs_metals_pc1.csv")

# Heatmap of top taxa vs contaminants (if enough)
sig_taxa <- gen_cor %>%
  dplyr::arrange(q_contam_pc1) %>%
  dplyr::filter(!is.na(q_contam_pc1)) %>%
  dplyr::slice_head(n = 40)

if (nrow(sig_taxa) >= 2) {
  mat <- sig_taxa %>%
    dplyr::select(rho_cu, rho_ni, rho_zn, rho_co, rho_contam_pc1) %>%
    as.matrix()
  rownames(mat) <- sig_taxa$Taxon
  png("exports/heatmap_taxa_rho_vs_metals_pc1.png", width=1400, height=1600, res=200)
  pheatmap::pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                     main = paste0("Taxa (", used_rank, ") ~ contaminants (Spearman rho)"))
  dev.off()
}

# ---------------------------
# Optional: ANCOM-BC2 (if installed & enough replication)
# ---------------------------
if (requireNamespace("ANCOMBC", quietly = TRUE)) {
  message("ANCOM-BC2: checking replication by contamination tier...")
  an_md <- as(phyloseq::sample_data(ps_coll), "data.frame")
  # Use contamination tier; require >=3 per group
  if ("contam_tier" %in% names(an_md)) {
    an_md$contam_tier <- factor(an_md$contam_tier, levels = c("Low","Medium","High"))
  } else {
    an_md$contam_tier <- NA
  }
  ok_rep <- all(table(an_md$contam_tier) >= 3, na.rm = TRUE) && nlevels(droplevels(an_md$contam_tier)) >= 2
  if (ok_rep) {
    out <- try(ANCOMBC::ancombc2(
      data = ps_coll,                # use collapsed counts (not relative)
      fix_formula = "contam_tier + zone",
      p_adj_method = "BH",
      prv_cut = 0.05, lib_cut = 1000,
      s0_perc = 0.05, alpha = 0.1, n_cl = 1
    ), silent = TRUE)
    if (!inherits(out, "try-error")) {
      tx_tbl <- as.data.frame(tax_table(ps_coll)) %>%
        tibble::rownames_to_column("Feature")
      da_tab <- out$res$diff_abn %>%
        tibble::rownames_to_column("Feature") %>%
        dplyr::left_join(tx_tbl, by = "Feature") %>%
        dplyr::arrange(q_val)
      readr::write_csv(da_tab, "exports/ancombc2_taxa_contam_tier.csv")
    } else {
      writeLines(paste0("ANCOM-BC2 skipped: ", as.character(out)),
                 "exports/ancombc2_taxa_contam_tier.txt")
    }
  } else {
    writeLines("ANCOM-BC2 skipped: need ≥3 samples per contamination tier (and ≥2 tiers present).",
               "exports/ancombc2_taxa_contam_tier.txt")
  }
} else {
  writeLines("ANCOM-BC2 not installed; skipping differential abundance.",
             "exports/ancombc2_taxa_contam_tier.txt")
}


# ---------------------------
# 9) Session info & README
# ---------------------------
readr::write_lines(c(
  "Vieques Microbiome Pipeline - Outputs overview",
  "exports/alpha_* : alpha diversity plots & stats",
  "exports/permanova_* : PERMANOVA results",
  "exports/betadisper_* : dispersion tests",
  "exports/ordination_* : ordination figures",
  "exports/capscale_* : dbRDA object",
  "exports/envfit_* : envfit loadings",
  "exports/genus_* : genus-level correlation/DA results",
  "exports/pcoa_scores_with_groups.csv : ordination scores table",
  "exports/join_diagnostics_samples_vs_chem.csv : sample↔chem join success table",
  "exports/chem_clean_with_index.csv : chemistry + contamination indices used"
), "exports/README_outputs.txt")

sink("exports/sessionInfo.txt"); print(sessionInfo()); sink()
message("✅ Pipeline complete. See the 'exports/' folder for results.")





# ===============================================================
# view_exports.R — quick viewer + HTML gallery for pipeline outputs
# ===============================================================

# Packages (only common, light deps)
need <- c("readr","dplyr","tibble","stringr","janitor","png","grid")
for(p in need){
  if(!requireNamespace(p, quietly = TRUE)){
    install.packages(p, repos="https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

exp_dir <- "exports"
stopifnot(dir.exists(exp_dir))

# Gather files
all_files <- list.files(exp_dir, full.names = TRUE)
pngs <- all_files[str_detect(tolower(all_files), "\\.png$")]
csvs <- all_files[str_detect(tolower(all_files), "\\.csv$")]
txts <- all_files[str_detect(tolower(all_files), "\\.txt$")]
rds  <- all_files[str_detect(tolower(all_files), "\\.rds$")]
others <- setdiff(all_files, c(pngs, csvs, txts, rds))

# ---------------------------
# 1) Console summary
# ---------------------------
cat("\n=== EXPORTS SUMMARY ===\n")
cat("PNG plots: ", length(pngs), "\n")
cat("CSV tables:", length(csvs), "\n")
cat("TXT files: ", length(txts), "\n")
cat("RDS files: ", length(rds), "\n")
cat("Other:     ", length(others), "\n\n")

# ---------------------------
# 2) Optional: show all PNGs in R plotting window
#     (comment this block if you don't want pop-ups)
# ---------------------------
if(length(pngs)){
  cat("Showing PNGs in the plotting window...\n")
  for (p in sort(pngs)) {
    try({
      img <- png::readPNG(p)
      grid::grid.newpage(); grid::grid.raster(img)
      title_txt <- basename(p)
      # draw a lightweight title
      grid::grid.text(title_txt, x=0.02, y=0.98, just=c("left","top"))
      Sys.sleep(0.5)
    }, silent = TRUE)
  }
}

# ---------------------------
# 3) Build an HTML gallery (exports/_index.html)
# ---------------------------
# helper: html escape
html_escape <- function(x){
  x <- gsub("&","&amp;",x, fixed=TRUE)
  x <- gsub("<","&lt;",x, fixed=TRUE)
  x <- gsub(">","&gt;",x, fixed=TRUE)
  x
}

# read first n lines of a text file
head_text <- function(path, n = 40){
  if(!file.exists(path)) return("")
  ln <- try(readLines(path, warn = FALSE), silent = TRUE)
  if(inherits(ln, "try-error")) return("")
  paste(utils::head(ln, n), collapse = "\n")
}

# read first rows of a CSV safely
head_csv_html <- function(path, n = 10){
  out <- try(suppressMessages(readr::read_csv(path, n_max = n, show_col_types = FALSE)), silent = TRUE)
  if(inherits(out, "try-error")) return("<em>(Could not preview CSV)</em>")
  out <- janitor::remove_empty(out, c("rows","cols"))
  # build simple HTML table
  if(nrow(out)==0) return("<em>(Empty CSV)</em>")
  hdr <- paste(sprintf("<th>%s</th>", html_escape(colnames(out))), collapse = "")
  rows <- apply(out, 1, function(r){
    cells <- paste(sprintf("<td>%s</td>", html_escape(as.character(r))), collapse = "")
    sprintf("<tr>%s</tr>", cells)
  })
  paste0("<table class='mini'><thead><tr>", hdr, "</tr></thead><tbody>",
         paste(rows, collapse = "\n"), "</tbody></table>")
}

# Minimal CSS
css <- "
body { font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif; margin: 24px; }
h1 { margin-top: 0; }
.grid { display: grid; gap: 16px; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); }
.card { border: 1px solid #e5e7eb; border-radius: 12px; padding: 12px; box-shadow: 0 1px 2px rgba(0,0,0,.03); }
.card img { max-width: 100%; height: auto; border-radius: 8px; display: block; }
.kv { color: #555; font-size: 0.9em; margin: 0 0 6px 0; }
a { color: #0b70ff; text-decoration: none; }
a:hover { text-decoration: underline; }
.mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; white-space: pre-wrap; background: #f8fafc; padding: 8px; border-radius: 8px; border: 1px solid #e5e7eb; }
.mini { border-collapse: collapse; width: 100%; font-size: 12px; }
.mini th, .mini td { border: 1px solid #e5e7eb; padding: 4px 6px; vertical-align: top; }
.section { margin-top: 28px; }
footer { margin-top: 24px; color: #666; font-size: 12px; }
"

html <- c(
  "<!doctype html>",
  "<html><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width, initial-scale=1'>",
  "<title>Exports Gallery</title>",
  "<style>", css, "</style>",
  "</head><body>",
  "<h1>Exports Gallery</h1>",
  sprintf("<p class='kv'>Directory: <span class='mono'>%s</span></p>", html_escape(normalizePath(exp_dir)))
)

# PNG section
html <- c(html, "<div class='section'><h2>Plots (PNG)</h2>")
if(length(pngs)){
  html <- c(html, "<div class='grid'>")
  for (p in sort(pngs)) {
    rel <- file.path(basename(exp_dir), basename(p))
    html <- c(html,
              "<div class='card'>",
              sprintf("<div class='kv'>%s</div>", html_escape(basename(p))),
              sprintf("<a href='%s'><img src='%s' alt='%s'></a>", html_escape(rel), html_escape(rel), html_escape(basename(p))),
              "</div>")
  }
  html <- c(html, "</div>")
} else {
  html <- c(html, "<p><em>No PNG plots found.</em></p>")
}
html <- c(html, "</div>")

# CSV section
html <- c(html, "<div class='section'><h2>Tables (CSV)</h2>")
if(length(csvs)){
  for (cpath in sort(csvs)) {
    rel <- file.path(basename(exp_dir), basename(cpath))
    html <- c(html,
              "<div class='card'>",
              sprintf("<div class='kv'><a href='%s'>%s</a></div>", html_escape(rel), html_escape(basename(cpath))),
              head_csv_html(cpath),
              "</div>")
  }
} else {
  html <- c(html, "<p><em>No CSV files found.</em></p>")
}
html <- c(html, "</div>")

# TXT section
html <- c(html, "<div class='section'><h2>Text Results (TXT)</h2>")
if(length(txts)){
  for (tpath in sort(txts)) {
    rel <- file.path(basename(exp_dir), basename(tpath))
    html <- c(html,
              "<div class='card'>",
              sprintf("<div class='kv'><a href='%s'>%s</a></div>", html_escape(rel), html_escape(basename(tpath))),
              sprintf("<div class='mono'>%s</div>", html_escape(head_text(tpath, 80))),
              "</div>")
  }
} else {
  html <- c(html, "<p><em>No TXT files found.</em></p>")
}
html <- c(html, "</div>")

# RDS section
html <- c(html, "<div class='section'><h2>Objects (RDS)</h2>")
if(length(rds)){
  html <- c(html, "<ul>")
  for (r in sort(rds)) {
    rel <- file.path(basename(exp_dir), basename(r))
    html <- c(html, sprintf("<li><a href='%s'>%s</a> — <span class='kv'>R object (not previewed)</span></li>",
                            html_escape(rel), html_escape(basename(r))))
  }
  html <- c(html, "</ul>")
} else {
  html <- c(html, "<p><em>No RDS files found.</em></p>")
}
html <- c(html, "</div>")

# Other files
if(length(others)){
  html <- c(html, "<div class='section'><h2>Other files</h2><ul>")
  for (o in sort(others)) {
    rel <- file.path(basename(exp_dir), basename(o))
    html <- c(html, sprintf("<li><a href='%s'>%s</a></li>", html_escape(rel), html_escape(basename(o))))
  }
  html <- c(html, "</ul></div>")
}

html <- c(html,
          "<footer>Generated by view_exports.R</footer>",
          "</body></html>"
)

# Write HTML
out_html <- file.path(exp_dir, "_index.html")
writeLines(html, out_html)
cat("✅ Wrote HTML gallery:", out_html, "\n")

# Open in default browser (comment out if undesired)
try(utils::browseURL(out_html), silent = TRUE)







# ===============================================================
# PCA of contaminants + Genus centroids overlay
# Produces: 
#  - exports/pca_contaminants_samples_w_genus_centroids.png
#  - exports/genus_centroids_on_contaminants_pca.csv
#  - exports/pca_contaminants_scores_by_sample.csv
#  - exports/pca_contaminants_loadings.csv
# ===============================================================

# ===============================================================
# PCA of contaminants + organics with OVERALL & BY-ZONE Genus centroids
# Outputs (saved in exports/):
#  - chemorg_pca_scores_by_site.csv
#  - chemorg_pca_loadings.csv
#  - chemorg_pca_samples.png
#  - chemorg_pca_samples_w_genus_centroids.png
#  - chemorg_pca_genus_centroids_overall.csv
#  - chemorg_pca_genus_centroids_by_zone.csv
#  - chemorg_pca_genus_centroids_by_zone.png
#  - chemorg_pca_biplot_loadings.png
# ===============================================================

# ---- deps# ===============================================================
# Flexible PCA plotting: by ZONE, by LOCATION, or SAMPLE ID
# Also: genus centroids overall / by zone / by location
# Files written under exports/
# ===============================================================

# --- helpers ---
.ensure_pkg <- function(p){ if(!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org"); suppressPackageStartupMessages(library(p, character.only=TRUE)) }
.ensure_pkg("ggrepel")

# ensures factor with >=2 levels to use as color; else returns NULL
.valid_group <- function(df, col){
  if (!col %in% names(df)) return(NULL)
  f <- droplevels(as.factor(df[[col]]))
  if (nlevels(f) >= 2) f else NULL
}

# Build a sample plot in different modes: "zone", "location", "samples"
plot_chem_pca <- function(sample_scores, pc1v, pc2v, mode = c("zone","location","samples"), outfile){
  mode <- match.arg(mode)
  df <- sample_scores
  has_zone <- !is.null(.valid_group(df, "zone"))
  has_loc  <- !is.null(.valid_group(df, "site_abbrev"))
  
  p <- ggplot(df, aes(PC1, PC2, label = SampleID))
  
  if (mode == "zone" && has_zone){
    p <- p + aes(color = zone, shape = site_abbrev)
    leg_sub <- c(color = "Zone", shape = "Location")
  } else if (mode == "location" && has_loc){
    p <- p + aes(color = site_abbrev)
    leg_sub <- c(color = "Location")
  } else {
    leg_sub <- NULL
  }
  
  p <- p +
    geom_point(size = 2.2, alpha = 0.9, na.rm = TRUE) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 60, show.legend = FALSE) +
    labs(
      title = paste("Contaminants + Organics PCA —", switch(mode, zone="by Zone", location="by Location", samples="per Sample")),
      x = glue::glue("PC1 ({pc1v}%)"), y = glue::glue("PC2 ({pc2v}%)")
    ) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  
  if (!is.null(leg_sub)){
    if ("color" %in% names(leg_sub) && "shape" %in% names(leg_sub)){
      p <- p + labs(color = leg_sub[["color"]], shape = leg_sub[["shape"]])
    } else if ("color" %in% names(leg_sub)){
      p <- p + labs(color = leg_sub[["color"]])
    }
  }
  
  ggsave(outfile, p, width = 9, height = 7, dpi = 300)
  p
}

# ------- SAMPLE-LEVEL PLOTS (three variants) -------
# Uses your existing sample_scores (from earlier block). 
# These will auto-fallback if zone/location is not informative.
plot_chem_pca(sample_scores, pc1v, pc2v, mode = "zone",     outfile = "exports/chemorg_pca_samples_by_zone.png")
plot_chem_pca(sample_scores, pc1v, pc2v, mode = "location", outfile = "exports/chemorg_pca_samples_by_location.png")
plot_chem_pca(sample_scores, pc1v, pc2v, mode = "samples",  outfile = "exports/chemorg_pca_samples_by_sampleid.png")

# ------- GENUS CENTROIDS: overall, by zone, by location -------
# long_abund built earlier: Feature, SampleID, rel_abund, Genus, site_abbrev, zone, PC1, PC2

# Overall centroid per Genus (weighted by rel_abund)
gen_centroids_overall <- long_abund |>
  dplyr::group_by(Genus) |>
  dplyr::summarise(
    w   = sum(rel_abund, na.rm = TRUE),
    PC1 = ifelse(w > 0, sum(PC1 * rel_abund, na.rm = TRUE) / w, NA_real_),
    PC2 = ifelse(w > 0, sum(PC2 * rel_abund, na.rm = TRUE) / w, NA_real_)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(is.finite(PC1), is.finite(PC2)) |>
  dplyr::arrange(dplyr::desc(w))
readr::write_csv(gen_centroids_overall, "exports/chemorg_pca_genus_centroids_overall.csv")

# By ZONE
gen_centroids_zone <- long_abund |>
  dplyr::group_by(Genus, zone) |>
  dplyr::summarise(
    w   = sum(rel_abund, na.rm = TRUE),
    PC1 = ifelse(w > 0, sum(PC1 * rel_abund, na.rm = TRUE) / w, NA_real_),
    PC2 = ifelse(w > 0, sum(PC2 * rel_abund, na.rm = TRUE) / w, NA_real_)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(is.finite(PC1), is.finite(PC2))
readr::write_csv(gen_centroids_zone, "exports/chemorg_pca_genus_centroids_by_zone.csv")

# By LOCATION
gen_centroids_loc <- long_abund |>
  dplyr::group_by(Genus, site_abbrev) |>
  dplyr::summarise(
    w   = sum(rel_abund, na.rm = TRUE),
    PC1 = ifelse(w > 0, sum(PC1 * rel_abund, na.rm = TRUE) / w, NA_real_),
    PC2 = ifelse(w > 0, sum(PC2 * rel_abund, na.rm = TRUE) / w, NA_real_)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(is.finite(PC1), is.finite(PC2))
readr::write_csv(gen_centroids_loc, "exports/chemorg_pca_genus_centroids_by_location.csv")




# ------- PLOTS with centroids (overall / by zone / by location) -------
topN <- 25
lab_overall <- gen_centroids_overall |> dplyr::slice_head(n = min(topN, n()))

# Overall
p_overall <- ggplot() +
  geom_point(data = sample_scores, aes(PC1, PC2, color = .valid_group(sample_scores, "zone") %||% NULL, shape = .valid_group(sample_scores, "site_abbrev") %||% NULL),
             size = 2, alpha = 0.85, na.rm = TRUE, show.legend = TRUE) +
  geom_point(data = gen_centroids_overall, aes(PC1, PC2, size = w), color = "black", alpha = 0.35, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = lab_overall, aes(PC1, PC2, label = Genus), size = 3, max.overlaps = 60, seed = 1) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  scale_size_continuous(name = "Overall genus abundance", range = c(1.2, 8)) +
  labs(title = glue::glue("Chem+Org PCA with Genus Centroids (top {nrow(lab_overall)} labeled)"),
       x = glue::glue("PC1 ({pc1v}%)"), y = glue::glue("PC2 ({pc2v}%)")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave("exports/chemorg_pca_samples_w_genus_centroids.png", p_overall, width = 9, height = 7, dpi = 300)

# By ZONE (facet)
keep_gen <- lab_overall$Genus
genz_plot <- gen_centroids_zone |> dplyr::filter(Genus %in% keep_gen)
p_by_zone <- ggplot() +
  geom_point(data = sample_scores, aes(PC1, PC2, color = zone), size = 1.8, alpha = 0.7, na.rm = TRUE) +
  geom_point(data = genz_plot, aes(PC1, PC2, size = w), color = "black", alpha = 0.45, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = genz_plot, aes(PC1, PC2, label = Genus), size = 2.8, seed = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  scale_size_continuous(name = "Genus abundance in zone", range = c(1.0, 7)) +
  labs(title = glue::glue("Genus centroids by Zone in Chem+Org PCA (top {length(keep_gen)} genera)"),
       x = glue::glue("PC1 ({pc1v}%)"), y = glue::glue("PC2 ({pc2v}%)")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  facet_wrap(~ zone, drop = TRUE)
ggsave("exports/chemorg_pca_genus_centroids_by_zone.png", p_by_zone, width = 11, height = 7.5, dpi = 300)

# By LOCATION (facet)
genl_plot <- gen_centroids_loc |> dplyr::filter(Genus %in% keep_gen)
p_by_loc <- ggplot() +
  geom_point(data = sample_scores, aes(PC1, PC2, color = site_abbrev), size = 1.8, alpha = 0.7, na.rm = TRUE) +
  geom_point(data = genl_plot, aes(PC1, PC2, size = w), color = "black", alpha = 0.45, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = genl_plot, aes(PC1, PC2, label = Genus), size = 2.5, seed = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  scale_size_continuous(name = "Genus abundance at location", range = c(1.0, 7)) +
  labs(title = glue::glue("Genus centroids by Location in Chem+Org PCA (top {length(keep_gen)} genera)"),
       x = glue::glue("PC1 ({pc1v}%)"), y = glue::glue("PC2 ({pc2v}%)")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  facet_wrap(~ site_abbrev, drop = TRUE)
ggsave("exports/chemorg_pca_genus_centroids_by_location.png", p_by_loc, width = 12, height = 8.5, dpi = 300)


message("✅ Chem+Org PCA with genus centroids complete. See files in exports/.")















######
######
######


# ===============================================================
# Phylogenetic tree + per-site heatmaps (top 30 taxa)
# ===============================================================

dir.create("exports", showWarnings = FALSE)
dir.create("exports/tree", showWarnings = FALSE)
dir.create("exports/heatmaps_by_site", showWarnings = FALSE)

# ---------- 1) Build/attach a phylogenetic tree ----------
# Use existing tree if present; otherwise create a taxonomy-based proxy tree.
if (is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))) {
  message("No phylogenetic tree detected; building a taxonomy-proxy tree.")
  tt <- as.data.frame(phyloseq::tax_table(ps))
  if (nrow(tt) < 2) stop("Need >=2 taxa to build a tree.")
  tt[is.na(tt)] <- ""
  ranks <- intersect(c("phylum","class","order","family","genus","species"), tolower(colnames(tt)))
  if (length(ranks) == 0) ranks <- colnames(tt)
  
  # distance = number of rank mismatches (simple, robust taxonomy proxy)
  n <- nrow(tt)
  D <- matrix(0, n, n, dimnames = list(rownames(tt), rownames(tt)))
  for (i in seq_len(n)) {
    for (j in i:n) {
      d <- sum(as.character(tt[i, ranks, drop=FALSE]) != as.character(tt[j, ranks, drop=FALSE]))
      D[i, j] <- D[j, i] <- d
    }
  }
  hc <- stats::hclust(as.dist(D), method = "average")
  tree_proxy <- ape::as.phylo(hc)
  phyloseq::phy_tree(ps) <- tree_proxy
}

# Save a quick tree plot of ALL taxa
png("exports/tree/tree_all_taxa.png", width = 1600, height = 1200, res = 200)
plot(phyloseq::phy_tree(ps), main = "Phylogenetic (proxy) tree — all taxa")
dev.off()

# ---------- 2) Prep: collapse to a readable taxonomic rank ----------
# Prefer Genus → Family → Order → Phylum → (ASV fallback)
tax_tbl <- as.data.frame(phyloseq::tax_table(ps))
rank_order <- c("Genus","genus","Family","family","Order","order","Phylum","phylum")
use_rank <- intersect(rank_order, colnames(tax_tbl))
use_rank <- if (length(use_rank)) use_rank[1] else NA_character_

ps_tx <- if (!is.na(use_rank)) suppressWarnings(phyloseq::tax_glom(ps, taxrank = use_rank, NArm = FALSE)) else ps
ps_tx_rel <- phyloseq::transform_sample_counts(ps_tx, function(x) if (sum(x) == 0) x else x / sum(x))

# Helper to label taxa nicely
label_for_taxa <- function(ps_obj, row_ids){
  tx <- as.data.frame(phyloseq::tax_table(ps_obj))
  if (!is.na(use_rank) && use_rank %in% names(tx)) {
    out <- as.character(tx[row_ids, use_rank])
    out[is.na(out) | !nzchar(out)] <- row_ids[is.na(out) | !nzchar(out)]
  } else {
    out <- row_ids
  }
  make.names(out, unique = TRUE)
}

# ---------- 3) Per-site heatmaps: top 30 taxa within each site ----------
sd <- as(phyloseq::sample_data(ps_tx_rel), "data.frame")
stopifnot("site_abbrev" %in% names(sd))
sites <- sort(unique(as.character(sd$site_abbrev)))

for (s in sites) {
  # subset to site
  ps_s <- phyloseq::subset_samples(ps_tx_rel, site_abbrev == s)
  ps_s <- phyloseq::prune_samples(phyloseq::sample_sums(ps_s) > 0, ps_s)
  if (phyloseq::nsamples(ps_s) < 1 || phyloseq::ntaxa(ps_s) < 2) {
    message("Skipping site ", s, " (not enough data).")
    next
  }
  
  # matrix (taxa x samples), pick top 30 by total rel. abundance IN THIS SITE
  M <- as(phyloseq::otu_table(ps_s), "matrix")
  if (!phyloseq::taxa_are_rows(ps_s)) M <- t(M)
  keep <- order(rowSums(M, na.rm = TRUE), decreasing = TRUE)[seq_len(min(30, nrow(M)))]
  M <- M[keep, , drop = FALSE]
  
  # row labels
  rownames(M) <- label_for_taxa(ps_s, rownames(M))
  
  # sample annotations (optional if these exist)
  ann <- as(phyloseq::sample_data(ps_s), "data.frame")
  ann$SampleID <- rownames(ann)
  ann_df <- ann %>%
    dplyr::transmute(
      SampleID,
      Zone = if ("zone" %in% names(ann)) ann$zone else NA,
      Organic = if ("organic_pct" %in% names(ann)) suppressWarnings(as.numeric(ann$organic_pct)) else NA_real_,
      Survival = if ("survival_pct" %in% names(ann)) suppressWarnings(as.numeric(ann$survival_pct)) else NA_real_
    ) %>% tibble::column_to_rownames("SampleID")
  
  # ensure column order aligned to annotation
  M <- M[, rownames(ann_df), drop = FALSE]
  
  # heatmap
  out_png <- file.path("exports/heatmaps_by_site", paste0("heatmap_top30_", s, ".png"))
  pheatmap::pheatmap(
    M,
    scale = "row",                        # emphasize within-taxon variation across samples
    clustering_distance_cols = "euclidean",
    clustering_method = "average",
    annotation_col = ann_df,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    main = paste0("Top 30 taxa — site ", s),
    filename = out_png,
    width = 11, height = 8
  )
}

# ---------- 4) (Optional) Tree of top 30 taxa per site ----------
# Saves a small tree per site for the same taxa shown in each heatmap.
for (s in sites) {
  ps_s <- phyloseq::subset_samples(ps_tx_rel, site_abbrev == s)
  if (phyloseq::nsamples(ps_s) < 1 || phyloseq::ntaxa(ps_s) < 2) next
  M <- as(phyloseq::otu_table(ps_s), "matrix"); if (!phyloseq::taxa_are_rows(ps_s)) M <- t(M)
  keep <- order(rowSums(M, na.rm = TRUE), decreasing = TRUE)[seq_len(min(30, nrow(M)))]
  taxa_keep <- rownames(M)[keep]
  
  # prune phyloseq to the top taxa for this site
  ps_top <- phyloseq::prune_taxa(taxa_keep, ps_s)
  
  # plot tree (proxy or real, depending on ps)
  out_tree <- file.path("exports/tree", paste0("tree_top30_", s, ".png"))
  png(out_tree, width = 1400, height = 1100, res = 200)
  plot(phyloseq::phy_tree(ps_top),
       main = paste0("Phylogenetic (proxy) tree — top 30 taxa at ", s))
  dev.off()
}

message("✅ Wrote: exports/tree/*.png and exports/heatmaps_by_site/heatmap_top30_*.png")



######
######
######
######
# ===============================================================
# Finish: tree plot + per-site heatmaps (top 30 taxa)
# (safe even if you've run parts already)
# ===============================================================

dir.create("exports", showWarnings = FALSE)
dir.create("exports/tree", showWarnings = FALSE)
dir.create("exports/heatmaps_by_site", showWarnings = FALSE)

# ---------- A) Make/attach a phylogenetic tree if missing ----------
if (is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))) {
  message("No phylogenetic tree detected; building a taxonomy-proxy tree.")
  tt <- as.data.frame(phyloseq::tax_table(ps))
  stopifnot(nrow(tt) >= 2)
  tt[is.na(tt)] <- ""
  # choose available ranks, case-insensitive
  ranks <- tolower(colnames(tt))
  ranks <- intersect(c("phylum","class","order","family","genus","species"), ranks)
  if (!length(ranks)) ranks <- tolower(colnames(tt))
  n <- nrow(tt)
  D <- matrix(0, n, n, dimnames = list(rownames(tt), rownames(tt)))
  for (i in seq_len(n)) {
    for (j in i:n) {
      d <- sum(as.character(tt[i, ranks, drop=FALSE]) != as.character(tt[j, ranks, drop=FALSE]))
      D[i, j] <- D[j, i] <- d
    }
  }
  hc <- stats::hclust(as.dist(D), method = "average")
  phyloseq::phy_tree(ps) <- ape::as.phylo(hc)
}

# Save a quick tree plot of ALL taxa (PNG)
png("exports/tree/tree_all_taxa.png", width = 1600, height = 1200, res = 200)
plot(phyloseq::phy_tree(ps), main = "Phylogenetic (proxy) tree — all taxa")
dev.off()

# ---------- B) Collapse to a readable rank & relative abundance ----------
tax_tbl <- as.data.frame(phyloseq::tax_table(ps))
rank_order <- c("Genus","genus","Family","family","Order","order","Phylum","phylum")
use_rank <- intersect(rank_order, colnames(tax_tbl))
use_rank <- if (length(use_rank)) use_rank[1] else NA_character_

ps_tx <- if (!is.na(use_rank)) suppressWarnings(phyloseq::tax_glom(ps, taxrank = use_rank, NArm = FALSE)) else ps
ps_tx_rel <- phyloseq::transform_sample_counts(ps_tx, function(x) if (sum(x) == 0) x else x / sum(x))

label_for_taxa <- function(ps_obj, row_ids){
  tx <- as.data.frame(phyloseq::tax_table(ps_obj))
  if (!is.na(use_rank) && use_rank %in% names(tx)) {
    out <- as.character(tx[row_ids, use_rank])
    out[is.na(out) | !nzchar(out)] <- row_ids[is.na(out) | !nzchar(out)]
  } else out <- row_ids
  make.names(out, unique = TRUE)
}

# ---------- C) Ensure site_abbrev exists in sample_data ----------
sd <- as(phyloseq::sample_data(ps_tx_rel), "data.frame")
if (!"site_abbrev" %in% names(sd)) {
  # try to infer from common alternatives
  guess <- intersect(names(sd), c("site","site_name","site_abbreviation","site_abbrev","site_abbr","abbrev","abbr"))
  stopifnot(length(guess) >= 1)
  sd$site_abbrev <- toupper(trimws(as.character(sd[[guess[1]]])))
  sample_data(ps_tx_rel)$site_abbrev <- sd$site_abbrev
}
sites <- sort(unique(as.character(sample_data(ps_tx_rel)$site_abbrev)))

# ---------- D) Per-site heatmaps: top 30 taxa ----------
for (s in sites) {
  ps_s <- phyloseq::subset_samples(ps_tx_rel, site_abbrev == s)
  ps_s <- phyloseq::prune_samples(phyloseq::sample_sums(ps_s) > 0, ps_s)
  if (phyloseq::nsamples(ps_s) < 1 || phyloseq::ntaxa(ps_s) < 2) {
    message("Skipping site ", s, " (not enough data)."); next
  }
  
  M <- as(phyloseq::otu_table(ps_s), "matrix")
  if (!phyloseq::taxa_are_rows(ps_s)) M <- t(M)
  
  keep <- order(rowSums(M, na.rm = TRUE), decreasing = TRUE)[seq_len(min(30, nrow(M)))]
  M <- M[keep, , drop = FALSE]
  rownames(M) <- label_for_taxa(ps_s, rownames(M))
  
  # Optional sample annotations if columns exist
  ann <- as(phyloseq::sample_data(ps_s), "data.frame")
  ann$SampleID <- rownames(ann)
  ann_df <- ann %>%
    dplyr::transmute(
      SampleID,
      Zone = if ("zone" %in% names(ann)) ann$zone else NA,
      Organic = if ("organic_pct" %in% names(ann)) suppressWarnings(as.numeric(ann$organic_pct)) else NA_real_,
      Survival = if ("survival_pct" %in% names(ann)) suppressWarnings(as.numeric(ann$survival_pct)) else NA_real_
    ) %>% tibble::column_to_rownames("SampleID")
  
  # align columns to annotation
  common_cols <- intersect(colnames(M), rownames(ann_df))
  M <- M[, common_cols, drop = FALSE]
  ann_df <- ann_df[common_cols, , drop = FALSE]
  
  out_png <- file.path("exports/heatmaps_by_site", paste0("heatmap_top30_", s, ".png"))
  pheatmap::pheatmap(
    M,
    scale = "row",
    clustering_distance_cols = "euclidean",
    clustering_method = "average",
    annotation_col = ann_df,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    main = paste0("Top 30 taxa — site ", s),
    filename = out_png,
    width = 11, height = 8
  )
}

# ---------- E) (Optional) per-site tree for those top 30 ----------
for (s in sites) {
  ps_s <- phyloseq::subset_samples(ps_tx_rel, site_abbrev == s)
  ps_s <- phyloseq::prune_samples(phyloseq::sample_sums(ps_s) > 0, ps_s)
  if (phyloseq::nsamples(ps_s) < 1 || phyloseq::ntaxa(ps_s) < 2) next
  
  M <- as(phyloseq::otu_table(ps_s), "matrix")
  if (!phyloseq::taxa_are_rows(ps_s)) M <- t(M)
  keep <- order(rowSums(M, na.rm = TRUE), decreasing = TRUE)[seq_len(min(30, nrow(M)))]
  ps_top <- phyloseq::prune_taxa(rownames(M)[keep], ps_s)
  
  out_tree <- file.path("exports/tree", paste0("tree_top30_", s, ".png"))
  png(out_tree, width = 1400, height = 1100, res = 200)
  plot(phyloseq::phy_tree(ps_top),
       main = paste0("Phylogenetic (proxy) tree — top 30 taxa at ", s))
  dev.off()
}

message("✅ Wrote: exports/tree/*.png and exports/heatmaps_by_site/heatmap_top30_*.png")

######
######
######
######

## ===================== PATCH: ORGANIC%, SURVIVAL%, ENV =====================

# Helper: robust numeric parser (handles numeric, factor, and "12.3%" strings)
to_num <- function(x){
  if (is.numeric(x)) return(as.numeric(x))
  if (is.factor(x))  x <- as.character(x)
  readr::parse_number(as.character(x))
}

# Helper: pick first existing column (define if not already present)
if (!exists("first_existing_col")) {
  first_existing_col <- function(df, cand){
    cc <- intersect(cand, names(df)); if (length(cc)) cc[1] else NULL
  }
}

# --- ORGANIC % from CHEM (e.g., "Percent Org Total") ---
# After clean_names(), likely becomes "percent_org_total"
ORGANIC_COL <- first_existing_col(
  chem,
  c("percent_org_total","percent_organic_total","percent_organic","percent_org","organic_percent")
)

# Create chem$organic_pct as numeric percent (3.21% -> 3.21)
chem$organic_pct <- if (!is.null(ORGANIC_COL)) to_num(chem[[ORGANIC_COL]]) else NA_real_

# Normalize key for the join
chem$site_abbrev <- toupper(trimws(as.character(chem$site_abbrev)))
mdf$site_abbrev  <- toupper(trimws(as.character(mdf$site_abbrev)))

# Merge organic % into MDF (by site_abbrev)
mdf <- dplyr::left_join(
  mdf,
  chem |>
    dplyr::select(site_abbrev, organic_pct,
                  cu = !!first_existing_col(chem, c("cu")), 
                  ni = !!first_existing_col(chem, c("ni")),
                  zn = !!first_existing_col(chem, c("zn")),
                  co = !!first_existing_col(chem, c("co"))),
  by = "site_abbrev"
)

# --- SURVIVAL % from MDF (e.g., "survival_100_percent_96hrs") ---
SURV_COL <- first_existing_col(
  mdf,
  c("survival_100_percent_96hrs","survival_percent","percent_survival","surv_pct","survival")
)

mdf$survival_pct <- if (!is.null(SURV_COL)) to_num(mdf[[SURV_COL]]) else NA_real_

# If survival is 0–1, convert to 0–100
if (is.finite(suppressWarnings(max(mdf$survival_pct, na.rm = TRUE))) &&
    suppressWarnings(max(mdf$survival_pct, na.rm = TRUE)) <= 1) {
  mdf$survival_pct <- mdf$survival_pct * 100
}

# --- Build env_keep safely (no rownames() inside dplyr) ---
mdf$SampleID <- if ("SampleID" %in% names(mdf) && !all(is.na(mdf$SampleID))) {
  as.character(mdf$SampleID)
} else {
  # fallback to rownames if needed
  if (!is.null(rownames(mdf))) rownames(mdf) else seq_len(nrow(mdf)) |> as.character()
}

env_keep <- mdf |>
  dplyr::transmute(
    SampleID,
    cu = suppressWarnings(as.numeric(cu)),
    ni = suppressWarnings(as.numeric(ni)),
    zn = suppressWarnings(as.numeric(zn)),
    co = suppressWarnings(as.numeric(co)),
    organic_pct = suppressWarnings(as.numeric(organic_pct)),
    survival_pct = suppressWarnings(as.numeric(survival_pct)),
    zone = droplevels(factor(zone)),
    site_abbrev = droplevels(factor(site_abbrev))
  )

# Diagnostics
message("✅ organic_pct non-NA: ", sum(is.finite(env_keep$organic_pct)))
message("✅ survival_pct non-NA: ", sum(is.finite(env_keep$survival_pct)))
readr::write_csv(env_keep, "exports/env_keep_for_CCA.csv")

# =================== END PATCH (continue with CCA/plots) ====================

# ===============================================================
# A) Detect and merge ORGANIC % and SURVIVAL % into metadata (mdf)
# ===============================================================

# --- Detect columns ---
ORGANIC_COL <- first_existing_col(chem, c("percent_org_total", "percent_organic", "percent_org", "percent_organics"))
SURV_COL    <- first_existing_col(mdf, c("survival_100_percent_96hrs", "survival_percent", "survival", "surv_96h"))

message("Detected organic column in chem: ", ORGANIC_COL %||% "None")
message("Detected survival column in sample data: ", SURV_COL %||% "None")

# --- Clean organic % from chem and merge into mdf ---
if(!is.null(ORGANIC_COL) && ORGANIC_COL %in% names(chem)){
  chem <- chem %>%
    mutate(
      !!ORGANIC_COL := readr::parse_number(!!sym(ORGANIC_COL))  # clean % signs to numeric
    )
  org_df <- chem %>% dplyr::select(site_abbrev, !!ORGANIC_COL)
  names(org_df)[2] <- "organic_pct"
  mdf <- mdf %>% dplyr::left_join(org_df, by = "site_abbrev")
} else {
  mdf$organic_pct <- NA_real_
}

# --- Clean survival % from sample data ---
mdf$survival_pct <- if(!is.null(SURV_COL) && SURV_COL %in% names(mdf)){
  readr::parse_number(mdf[[SURV_COL]])
} else NA_real_

# ===============================================================
# B) Build safe environmental matrix for CCA and ordinations
# ===============================================================

env_keep <- mdf %>%
  dplyr::transmute(
    SampleID = rownames(.),
    cu = as.numeric(cu),
    ni = as.numeric(ni),
    zn = as.numeric(zn),
    co = as.numeric(co),
    organic_pct = as.numeric(organic_pct),
    survival_pct = as.numeric(survival_pct),
    zone = droplevels(factor(zone)),
    site_abbrev = droplevels(factor(site_abbrev))
  )

# ===============================================================
# C) Confirm completeness
# ===============================================================
summary(env_keep)
readr::write_csv(env_keep, "exports/env_keep_for_CCA.csv")
message("✅ Environmental metadata ready for CCA and vector-based analyses.")

# ===============================================================
# D) Canonical Correspondence Analysis (CCA)
# ===============================================================

# Select taxa table (relative abundance)
ps_rel <- phyloseq::transform_sample_counts(ps, function(x) if(sum(x)==0) x else x/sum(x))
otu_rel <- as.data.frame(phyloseq::otu_table(ps_rel))
otu_rel <- otu_rel[, match(env_keep$SampleID, colnames(otu_rel)), drop = FALSE]

# Run CCA using contaminants + organics + survival
cca_model <- vegan::cca(t(otu_rel) ~ cu + ni + zn + co + organic_pct + survival_pct, data = env_keep)

# Plot CCA biplot
png("exports/cca_contaminants_organics_survival.png", width=1200, height=900, res=200)
plot(cca_model, display=c("sites","species","bp"), main="CCA — Contaminants, Organics & Survival")
dev.off()

# Export scores
scores_sites <- as.data.frame(scores(cca_model, display="sites"))
scores_species <- as.data.frame(scores(cca_model, display="species"))
readr::write_csv(scores_sites, "exports/cca_scores_sites.csv")
readr::write_csv(scores_species, "exports/cca_scores_species.csv")

# ===============================================================
# E) Phylogenetic Tree by Sample
# ===============================================================

tree <- phyloseq::phy_tree(ps)
if(!is.null(tree)){
  png("exports/phylogenetic_tree_by_sample.png", width=1400, height=1000, res=200)
  plot(tree, main="Phylogenetic Tree of Taxa per Sample")
  dev.off()
} else {
  message("⚠️ No phylogenetic tree in phyloseq object.")
}

# ===============================================================
# F) Heatmap: Survival × Contaminants × Dominant Taxa
# ===============================================================

# Identify most prevalent taxa
taxa_abund <- data.frame(taxa = taxa_names(ps), abundance = taxa_sums(ps))
top_taxa <- taxa_abund %>% arrange(desc(abundance)) %>% slice_head(n = 40) %>% pull(taxa)

# Create matrix of top taxa (relative abundance)
ps_top <- prune_taxa(top_taxa, ps_rel)
mat <- as.data.frame(otu_table(ps_top))
mat <- mat %>% mutate(Taxon = rownames(.)) %>% relocate(Taxon)

# Combine with survival and contaminants
meta_for_heat <- env_keep %>%
  dplyr::select(SampleID, cu, ni, zn, co, organic_pct, survival_pct)
rownames(meta_for_heat) <- meta_for_heat$SampleID

# Order samples by survival
mat_order <- mat[, c("Taxon", meta_for_heat$SampleID[order(meta_for_heat$survival_pct, decreasing = TRUE)])]

pheatmap::pheatmap(
  mat_order %>% column_to_rownames("Taxon"),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  annotation_col = meta_for_heat %>% select(survival_pct, organic_pct),
  main = "Heatmap: Top Taxa vs Contaminants & Survival"
)

# ===============================================================
# G) Alpha & Beta Diversity by Location
# ===============================================================

# Alpha diversity boxplots
alpha_df <- readr::read_csv("exports/alpha_metrics_merged.csv", show_col_types = FALSE)
p_alpha_loc <- ggplot(alpha_df, aes(x = site_abbrev, y = Shannon, fill = site_abbrev)) +
  geom_boxplot() + geom_jitter(width = 0.2) +
  labs(title = "Alpha Diversity (Shannon) by Location", x = "Location", y = "Shannon Index") +
  theme_bw(base_size = 12) + theme(legend.position="none")
ggsave("exports/alpha_diversity_by_location.png", p_alpha_loc, width=8, height=5, dpi=300)

# Beta diversity PCoA colored by SampleID
ord_pcoa <- ordinate(ps_rel, method="PCoA", distance="bray")
scores <- as.data.frame(ord_pcoa$vectors) %>% tibble::rownames_to_column("SampleID")
p_beta <- ggplot(scores, aes(Axis.1, Axis.2, color = SampleID)) +
  geom_point(size=2.2, alpha=0.8) +
  ggrepel::geom_text_repel(aes(label=SampleID), size=2.5, show.legend=FALSE) +
  theme_bw(base_size=12) + theme(legend.position="none") +
  labs(title="Beta Diversity (PCoA — Bray-Curtis)", x="PCoA1", y="PCoA2")
ggsave("exports/beta_diversity_by_sampleid.png", p_beta, width=8, height=6, dpi=300)

# ===============================================================
# ✅ END — Ready for integration into your main pipeline
# ===============================================================



######
######

# ---------- A) Identify organic % and survival columns ----------
# From chemistry: try to find an "organic %" proxy (e.g., LOI/TOC/etc.)
find_org_col <- function(df){
  cand <- names(df)
  # prioritize % or fractional organics names
  keys <- c("organic", "org", "loi", "toc", "t.o.c", "total_organic", "percent_organic", "%_organic")
  hits <- cand[Reduce(`|`, lapply(keys, function(k) grepl(k, cand, ignore.case = TRUE)))]
  # avoid false positives like "organization"
  hits <- hits[!grepl("organizat|organis", hits, ignore.case = TRUE)]
  hits[1]
}

# From metadata: try to find survival %
find_surv_col <- function(df){
  cand <- names(df)
  keys <- c("survival", "surv", "mortality", "percent_survival", "surv_pct", "%survival")
  hits <- cand[Reduce(`|`, lapply(keys, function(k) grepl(k, cand, ignore.case = TRUE)))]
  hits[1]
}

ORGANIC_COL <- find_org_col(chem)
SURV_COL    <- find_surv_col(mdf)

message("Detected organic column in chem: ", ORGANIC_COL %||% "None")
message("Detected survival column in sample data: ", SURV_COL %||% "None")

# Coerce numeric copies on mdf for modeling/annotations
mdf$organic_pct <- if (!is.null(ORGANIC_COL) && ORGANIC_COL %in% names(chem)) {
  # merge chem organic into mdf by site_abbrev
  org_df <- chem |> dplyr::select(site_abbrev, !!ORGANIC_COL)
  names(org_df)[2] <- "organic_raw"
  tmp <- mdf |> dplyr::left_join(org_df, by = "site_abbrev")
  suppressWarnings(as.numeric(tmp$organic_raw))
} else NA_real_

mdf$survival_pct <- if (!is.null(SURV_COL) && SURV_COL %in% names(mdf)) {
  suppressWarnings(as.numeric(mdf[[SURV_COL]]))
} else NA_real_

# A safe “env” frame we’ll reuse
env_keep <- mdf |>
  dplyr::transmute(
    SampleID = rownames(.),
    cu = as.numeric(cu), ni = as.numeric(ni), zn = as.numeric(zn), co = as.numeric(co),
    organic_pct = as.numeric(organic_pct),
    survival_pct = as.numeric(survival_pct),
    zone = droplevels(factor(zone)),
    site_abbrev = droplevels(factor(site_abbrev))
  )




# ====
# 

# ---------- B) CCA with arrows & labeled samples ----------
# Hellinger-transform the community (recommended for CCA)
ps_rel0 <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
comm    <- as(otu_table(ps_rel0), "matrix"); if(!taxa_are_rows(ps_rel0)) comm <- t(comm)
comm_hel <- vegan::decostand(comm, method = "hellinger")

# Align env to samples in comm
env2 <- env_keep |> dplyr::filter(SampleID %in% colnames(comm_hel))
comm_hel <- comm_hel[, env2$SampleID, drop = FALSE]

# Select numeric predictors that exist
vars_num <- c("cu","ni","zn","co","organic_pct","survival_pct")
X <- env2[, vars_num, drop = FALSE] |> dplyr::mutate(dplyr::across(everything(), as.numeric))

# Remove columns with all NA or zero variance
good_cols <- names(X)[sapply(X, function(v) sum(is.finite(v)) >= 3 && stats::var(v, na.rm=TRUE) > 0)]
X <- X[, good_cols, drop = FALSE]

stopifnot(nrow(env2) >= 3, ncol(comm_hel) >= 3, ncol(X) >= 1)

cca_fit <- vegan::cca(t(comm_hel) ~ ., data = X)  # sites = samples
# Save a quick summary
sink("exports/cca_summary.txt"); print(summary(cca_fit)); sink()

# Tidy site scores for plotting & labels
sc_sites <- as.data.frame(vegan::scores(cca_fit, display = "sites"))
sc_sites$SampleID <- rownames(sc_sites)
sc_sites <- sc_sites |>
  dplyr::left_join(mdf |> dplyr::select(SampleID = rownames(mdf), zone, site_abbrev), by = "SampleID")

# Tidy biplot arrows (env loadings)
sc_arrows <- as.data.frame(vegan::scores(cca_fit, display = "bp"))
sc_arrows$var <- rownames(sc_arrows)

library(ggrepel)
p_cca <- ggplot(sc_sites, aes(CCA1, CCA2, color = zone)) +
  geom_point(size = 2.4, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = 80, show.legend = FALSE) +
  geom_segment(data = sc_arrows,
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.6) +
  geom_text(data = sc_arrows, aes(CCA1, CCA2, label = var),
            vjust = -0.4, size = 3) +
  labs(title = "CCA — community vs contaminants, organic %, survival",
       color = "Zone") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave("exports/CCA_biplot_samples_labeled.png", p_cca, width = 9, height = 7, dpi = 300)

# 
# 
# 
## ---------- C) Tree of taxa with per-sample abundance heat (one panel per sample) ----------
# If a true phylogeny exists in `ps`, use it; else build a taxonomy-based proxy tree.
ps_tree <- ps
if (is.null(phyloseq::phy_tree(ps_tree, errorIfNULL = FALSE))) {
  message("No phylogenetic tree present; constructing a taxonomy-based proxy tree.")
  tt <- as.data.frame(phyloseq::tax_table(ps_tree))
  # Build a simple hierarchical string and cluster
  tax_path <- apply(tt, 1, function(x) paste(na.omit(as.character(x)), collapse = ";"))
  D <- stats::dist(stringdist::stringdistmatrix(tax_path, tax_path, method = "lv"))
  hc <- stats::hclust(D, method = "average")
  tree_proxy <- ape::as.phylo(hc)
  phyloseq::phy_tree(ps_tree) <- tree_proxy
}

# Make relative abundances and cap small values
ps_tree_rel <- phyloseq::transform_sample_counts(ps_tree, function(x) x/sum(x))
otu_rel <- as(otu_table(ps_tree_rel), "matrix"); if(!taxa_are_rows(ps_tree_rel)) otu_rel <- t(otu_rel)
otu_rel <- pmin(otu_rel, 0.2)

# Plot a multi-sample tree heat (one big heat), then a faceted PDF for per-sample
dir.create("exports/tree", showWarnings = FALSE)

# Overall tree with a compact heatmap strip per sample
pheatmap::pheatmap(
  otu_rel,
  show_rownames = FALSE, show_colnames = TRUE,
  filename = "exports/tree/taxa_tree_heat_overall.png", width = 14, height = 10
)

# Optional: one PNG per sample highlighting that sample’s abundance
samples_list <- colnames(otu_rel)
for (s in samples_list) {
  mat <- cbind(Abundance = otu_rel[, s, drop = FALSE])
  pheatmap::pheatmap(
    mat, show_rownames = FALSE, show_colnames = TRUE,
    filename = file.path("exports/tree", paste0("taxa_tree_heat_", s, ".png")),
    width = 5, height = 8
  )
}



###
###
###
# ---------- D) Heatmap of survival, contaminants, and top taxa ----------
# Get relative abundance by Genus (collapse safely)
rk <- intersect(c("Genus","genus","Family","family"), colnames(as.data.frame(tax_table(ps))))
use_rank <- if (length(rk)) rk[1] else NA_character_
ps_gen <- if (!is.na(use_rank)) suppressWarnings(tax_glom(ps, taxrank = use_rank, NArm = FALSE)) else ps
ps_gen_rel <- transform_sample_counts(ps_gen, function(x) x/sum(x))
otu_g <- as.data.frame(otu_table(ps_gen_rel)); if(!taxa_are_rows(ps_gen_rel)) otu_g <- t(otu_g)

# Pick top taxa by overall prevalence (sum of rel. abund)
topN <- 25
keep <- order(rowSums(otu_g, na.rm = TRUE), decreasing = TRUE)[seq_len(min(topN, nrow(otu_g)))]
mat_top <- as.matrix(otu_g[keep, , drop = FALSE])

# Taxon labels
tx <- as.data.frame(tax_table(ps_gen_rel))
lab <- if (!is.na(use_rank) && use_rank %in% names(tx)) tx[rownames(mat_top), use_rank] else rownames(mat_top)
rownames(mat_top) <- make.names(lab, unique = TRUE)

# Annotations for samples
ann <- mdf |>
  dplyr::transmute(
    SampleID = rownames(.),
    Survival = survival_pct,
    Cu = as.numeric(cu), Ni = as.numeric(ni), Zn = as.numeric(zn), Co = as.numeric(co),
    Organic = as.numeric(organic_pct),
    Zone = zone, Location = site_abbrev
  ) |>
  tibble::column_to_rownames("SampleID")

mat_top <- mat_top[, rownames(ann), drop = FALSE]  # ensure same order

pheatmap::pheatmap(
  mat_top,
  annotation_col = ann[, c("Survival","Cu","Ni","Zn","Co","Organic","Zone","Location")],
  scale = "row",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  show_rownames = TRUE,
  filename = "exports/heatmap_top_taxa_survival_contaminants.png",
  width = 14, height = 9
)


# 

# ---------- E) Alpha by location (site_abbrev) ----------
alpha_df <- readr::read_csv("exports/alpha_metrics_merged.csv", show_col_types = FALSE)
alpha_df$site_abbrev <- factor(alpha_df$site_abbrev)

p_alpha_loc <- ggplot(alpha_df, aes(site_abbrev, Shannon, fill = site_abbrev)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  labs(x = "Location (site_abbrev)", y = "Shannon",
       title = "Alpha diversity (Shannon) by location") +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  coord_flip()
ggsave("exports/alpha_Shannon_by_location.png", p_alpha_loc, width = 7, height = 6, dpi = 300)

# 
# ---------- F) PCoA with sample labels ----------
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))
dist_bray <- phyloseq::distance(ps_rel, method = "bray")
ord      <- ordinate(ps_rel, method = "PCoA", distance = dist_bray)

scores <- as.data.frame(ord$vectors)
scores$SampleID <- rownames(scores)
scores <- scores |>
  dplyr::left_join(mdf |> dplyr::select(SampleID = rownames(mdf), zone, site_abbrev), by = "SampleID")

pc1 <- round(ord$values$Relative_eig[1]*100,1)
pc2 <- round(ord$values$Relative_eig[2]*100,1)

p_lab <- ggplot(scores, aes(Axis.1, Axis.2, color = zone, label = SampleID)) +
  geom_point(size = 2.4, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 80, show.legend = FALSE) +
  labs(x = glue::glue("PCoA1 ({pc1}%)"), y = glue::glue("PCoA2 ({pc2}%)"),
       title = "PCoA (Bray-Curtis) — labeled by sample ID") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave("exports/beta_pcoa_labeled_by_sample.png", p_lab, width = 9, height = 7, dpi = 300)



###
###
###
# ---------- G) PCA on top taxa with envfit arrows ----------
# Use the same top taxa matrix from D (mat_top), apply Hellinger transform again just in case
mat_pca <- vegan::decostand(t(mat_top), method = "hellinger")  # samples x taxa
pca_fit <- stats::prcomp(mat_pca, center = TRUE, scale. = FALSE)

# Site scores
pca_sites <- as.data.frame(pca_fit$x) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::left_join(env_keep, by = "SampleID")

# Variance %
varp <- round(100 * summary(pca_fit)$importance["Proportion of Variance", 1:2], 1)

# Fit env arrows (contaminants + organic)
env_for_fit <- pca_sites |> dplyr::select(cu,ni,zn,co,organic_pct)
ef <- vegan::envfit(as.matrix(pca_fit$x[,1:2]), env_for_fit, permutations = 9999)

ef_ar <- as.data.frame(ef$vectors$arrows * rep(ef$vectors$r, each = 2))
names(ef_ar) <- c("PC1","PC2")
ef_ar$var <- rownames(ef$vectors$arrows)

p_taxa_pca <- ggplot(pca_sites, aes(PC1, PC2, color = zone)) +
  geom_point(size = 2.4, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = 80, show.legend = FALSE) +
  geom_segment(data = ef_ar,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.6) +
  geom_text(data = ef_ar, aes(PC1, PC2, label = var), vjust = -0.4, size = 3) +
  labs(title = "PCA of most prevalent taxa (+ envfit: contaminants & organic)",
       x = glue::glue("PC1 ({varp[1]}%)"), y = glue::glue("PC2 ({varp[2]}%)"),
       color = "Zone") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave("exports/pca_top_taxa_envfit.png", p_taxa_pca, width = 9, height = 7, dpi = 300)

# ===========================================================


# ===============================================================
# Vieques Microbiome Pipeline v4 — ASV table + taxonomy + metadata
# Author: auto-generated for Valeria
# Inputs (edit the paths if needed):
#   - "asv_table.csv"
#   - "asv_id_map_with_taxonomy.csv"                (ID column = 'asv')
#   - "Meta data  - Copy of sediment toxicity .csv" (contains SampleID)
#   - "CHEM_SED.xlsx"
# Outputs: exports/ folder with CSV tables, plots (PNG), and diagnostics.
# ===============================================================

# ---------------------------
# 0) Packages & helpers
# ---------------------------
need <- c(
  "tidyverse","data.table","janitor","readxl","stringr","glue",
  "vegan","phyloseq","ape","matrixStats","broom",
  "ggpubr","FSA","rstatix","pheatmap","scales","magrittr","tibble"
)
for(p in need){
  if(!requireNamespace(p, quietly = TRUE)){
    install.packages(p, repos="https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}
set.seed(42)
theme_set(theme_bw(base_size = 12))
dir.create("exports", showWarnings = FALSE)

# Helpers
stdz <- function(x) as.numeric(scale(x))
lz   <- function(x) stdz(log10(x + 1))
first_existing_col <- function(df, candidates) {
  cand <- intersect(candidates, names(df))
  if (length(cand)) cand[1] else NULL
}
coalesce_across <- function(df, candidates) {
  cand <- intersect(candidates, names(df))
  if (!length(cand)) return(rep(NA_character_, nrow(df)))
  out <- df[[cand[1]]]
  if (length(cand) > 1) {
    for (i in cand[-1]) {
      out <- dplyr::coalesce(dplyr::na_if(out, ""), dplyr::na_if(df[[i]], ""))
    }
  }
  out
}
clean_id <- function(x){
  x %>% as.character() %>% trimws() %>% str_replace_all("\\s+", "_")
}

# ---------------------------
# 1) User paths (EDIT if needed)
# ---------------------------
asv_path  <- "asv_table.csv"
tax_path  <- "asv_id_map_with_taxonomy.csv"
meta_path <- "Meta data  - Copy of sediment toxicity .csv"  # spaces intentional
chem_path <- "CHEM_SED.xlsx"

stopifnot(file.exists(asv_path), file.exists(tax_path),
          file.exists(meta_path), file.exists(chem_path))
message("✅ Inputs found.")

# ---------------------------
# 2) Chemistry → indices (PC1, z-mean, tiers)
# ---------------------------
chem_raw <- readxl::read_excel(chem_path) |> janitor::clean_names()

chem <- chem_raw |>
  rename(
    site = any_of(c("site","site_name")),
    zone = any_of(c("zone")),
    site_abbrev = any_of(c("site_abbreviation","site_abbrev","site_abbr","site__abbreviation",
                           "site__abbr","site_abbreviation_")),
    cu = any_of(c("cu_mg_kg","cu_mgkg","cu","cu_mg_k_g")),
    ni = any_of(c("ni_mg_kg","ni_mgkg","ni")),
    zn = any_of(c("zn_mg_kg","zn_mgkg","zn")),
    co = any_of(c("co_mg_kg","co_mgkg","co"))
  )

# fallback metal grabs if weird names
nm <- names(chem_raw)
grab_metal <- function(keys){
  cand <- nm[str_detect(nm, keys)]
  if(length(cand)) chem_raw[[cand[1]]] else NA_real_
}
if(!"cu" %in% names(chem)) chem$cu <- grab_metal("cu.*mg.*kg|\\bcu\\b")
if(!"ni" %in% names(chem)) chem$ni <- grab_metal("ni.*mg.*kg|\\bni\\b")
if(!"zn" %in% names(chem)) chem$zn <- grab_metal("zn.*mg.*kg|\\bzn\\b")
if(!"co" %in% names(chem)) chem$co <- grab_metal("co.*mg.*kg|\\bco\\b")

# site_abbrev derive if needed
if(is.null(chem$site_abbrev)){
  chem$site_abbrev <- coalesce_across(chem_raw,
                                      c("site_abbreviation","site_abbrev","abbrev","abbr","site_abbr","site__abbreviation"))
  if(is.null(chem$site_abbrev) && "site" %in% names(chem)) {
    chem$site_abbrev <- str_replace_all(toupper(substr(chem$site,1,4)), "[^A-Z0-9]", "")
  }
}

chem <- chem |>
  mutate(across(c(cu,ni,zn,co), as.numeric)) |>
  select(site, zone, site_abbrev, cu, ni, zn, co) |>
  distinct()

metals_trans <- chem |>
  select(cu,ni,zn,co) |>
  mutate(across(everything(), lz))

chem$contam_index_zmean <- rowMeans(metals_trans, na.rm = TRUE)
pc <- prcomp(metals_trans |> mutate(across(everything(), ~replace_na(.,0))),
             center = TRUE, scale. = FALSE)
chem$contam_pc1 <- pc$x[,1]

q <- quantile(chem$contam_pc1, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE) |> unique()
if (length(q) < 4) { rng <- range(chem$contam_pc1, na.rm = TRUE); q <- seq(rng[1], rng[2], length.out = 4) }
chem$contam_tier <- cut(chem$contam_pc1, breaks = q, include.lowest = TRUE,
                        labels = c("Low","Medium","High"))

readr::write_csv(chem, "exports/chem_clean_with_index.csv")

# ---------------------------
# 3) Metadata (load first; normalize IDs)
# ---------------------------
meta <- suppressMessages(readr::read_csv(meta_path, guess_max = 100000)) %>% janitor::clean_names()
sample_col <- first_existing_col(meta, c("sampleid","sample_id","sample","id","sample_name"))
if (is.null(sample_col)) stop("SampleID column not found in metadata. Available: ", paste(names(meta), collapse=", "))
names(meta)[names(meta)==sample_col] <- "SampleID"
meta$SampleID <- clean_id(meta$SampleID)

site_candidates <- c("site_abbrev","site_abbreviation","site_abbr","abbrev","abbr","site","site_name")
zone_col <- first_existing_col(meta, c("zone","Zone"))

meta <- meta %>%
  mutate(
    site_abbrev = coalesce_across(., site_candidates),
    site_abbrev = toupper(trimws(as.character(site_abbrev))),
    zone = if (!is.null(zone_col)) .[[zone_col]] else NA_character_
  )

# ---------------------------
# 4) ASV table + Taxonomy → phyloseq
# ---------------------------
# ASV table: auto-detect orientation and transpose if samples are rows
asv_raw <- suppressMessages(readr::read_csv(asv_path))
first_col <- names(asv_raw)[1]
looks_like_taxa   <- grepl("^(ASV|Feature|OTU|Seq|Sequence)", first_col, ignore.case = TRUE)
looks_like_sample <- grepl("^(Sample|SampleID|ID|Run|Library)", first_col, ignore.case = TRUE)

if (looks_like_sample && !looks_like_taxa) {
  message("Detected samples in rows; transposing ASV table so taxa are rows.")
  mat <- as.matrix(asv_raw[,-1])
  rownames(mat) <- clean_id(asv_raw[[1]])
  mat[is.na(mat)] <- 0
  mat <- t(mat)
  asv_df <- as.data.frame(mat) %>% tibble::rownames_to_column("ASV_ID")
} else {
  asv_df <- asv_raw
  if (!grepl("ASV|Feature|OTU|Seq", first_col, ignore.case = TRUE)) names(asv_df)[1] <- "ASV_ID"
  asv_df[[1]] <- clean_id(asv_df[[1]])
}
# numeric counts
for (j in 2:ncol(asv_df)) asv_df[[j]] <- suppressWarnings(as.numeric(asv_df[[j]]))
asv_df[is.na(asv_df)] <- 0

# Taxonomy: your file key column is 'asv'
tax_raw <- suppressMessages(readr::read_csv(tax_path)) %>% janitor::clean_names()
if (!"asv" %in% names(tax_raw)) stop("Taxonomy file must contain column 'asv'. Found: ", paste(names(tax_raw), collapse=", "))
tax_raw$asv <- clean_id(tax_raw$asv)

tax_df <- tax_raw %>%
  rename(ASV_ID = asv) %>%
  select(ASV_ID, any_of(c("kingdom","phylum","class","order","family","genus","species"))) %>%
  distinct()

# Build matrices
asv_ids <- asv_df$ASV_ID
otu_mat <- as.matrix(asv_df[,-1]); rownames(otu_mat) <- asv_ids

tax_mat <- tax_df %>%
  right_join(tibble(ASV_ID = asv_ids), by = "ASV_ID") %>%
  select(-ASV_ID) %>% as.matrix()
rownames(tax_mat) <- asv_ids

# Harmonize samples with metadata
colnames(otu_mat) <- clean_id(colnames(otu_mat))
meta$SampleID     <- clean_id(meta$SampleID)
common_samples <- intersect(colnames(otu_mat), meta$SampleID)
message(glue::glue("Samples in ASV table:   {ncol(otu_mat)}"))
message(glue::glue("Samples in metadata:    {nrow(meta)}"))
message(glue::glue("Overlap (by SampleID):  {length(common_samples)}"))
if (length(common_samples) < 2) stop("Fewer than 2 overlapping samples between ASV table and metadata after ID cleaning.")

# Subset & order
otu_mat <- otu_mat[, common_samples, drop=FALSE]
meta2   <- meta %>% filter(SampleID %in% common_samples) %>% distinct(SampleID, .keep_all = TRUE) %>% as.data.frame()
rownames(meta2) <- meta2$SampleID
meta2$zone <- as.factor(meta2$zone)

# phyloseq object
OTU  <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX  <- tax_table(tax_mat)
META <- sample_data(meta2)
ps   <- phyloseq(OTU, TAX, META)

# prune zeros
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)
stopifnot(nsamples(ps) > 1, ntaxa(ps) > 1)

# ---------------------------
# 5) Join chemistry onto samples (by site_abbrev)
# ---------------------------
# ---------------------------
# 5) Join chemistry onto samples (by site_abbrev)
# ---------------------------

# 5a) Start from a plain data.frame, not sample_data()
sd_df <- as(sample_data(ps), "data.frame")
sd_df$SampleID <- rownames(sd_df)

# Standardize key in chemistry
chem_merge <- chem %>%
  mutate(site_abbrev = toupper(trimws(as.character(site_abbrev))))

# 5b) Left-join chem → samples by site_abbrev
sd2_df <- sd_df %>%
  dplyr::left_join(chem_merge, by = "site_abbrev")

# Report join success
message(glue::glue(
  "Joined chemistry for {sum(!is.na(sd2_df$contam_pc1))}/{nrow(sd2_df)} samples (by site_abbrev)."
))

# ---- 5c’) Robustly synthesize a single `zone` column, then write CSV and assign ----

# Start with a blank character vector of the right length
sd2_df$zone <- rep(NA_character_, nrow(sd2_df))

# If zone.x exists (from metadata), prefer it when non-empty
if ("zone.x" %in% names(sd2_df)) {
  z1 <- as.character(sd2_df$zone.x)
  keep1 <- !is.na(z1) & nzchar(z1)
  sd2_df$zone[keep1] <- z1[keep1]
}

# If zone.y exists (from chemistry), fill any remaining NAs/empties
if ("zone.y" %in% names(sd2_df)) {
  z2 <- as.character(sd2_df$zone.y)   # coerce to character to avoid type clashes
  need2 <- is.na(sd2_df$zone) | !nzchar(sd2_df$zone)
  sd2_df$zone[need2] <- z2[need2]
}

# Make zone a factor for downstream models
sd2_df$zone <- factor(sd2_df$zone)

# Tidy up: drop the suffix columns if present
sd2_df <- sd2_df %>% dplyr::select(-dplyr::any_of(c("zone.x","zone.y")))

# Diagnostics CSV
readr::write_csv(
  sd2_df %>% dplyr::transmute(SampleID, site_abbrev, zone, has_chem = !is.na(contam_pc1)),
  "exports/join_diagnostics_samples_vs_chem.csv"
)

# Assign back to phyloseq
sd2_ps <- sd2_df
rownames(sd2_ps) <- sd2_ps$SampleID
sample_data(ps) <- phyloseq::sample_data(sd2_ps)

# ---- 5c’) Robustly synthesize a single `zone` column, then write CSV and assign ----

# Start with a blank character vector of the right length
sd2_df$zone <- rep(NA_character_, nrow(sd2_df))

# If zone.x exists (from metadata), prefer it when non-empty
if ("zone.x" %in% names(sd2_df)) {
  z1 <- as.character(sd2_df$zone.x)
  keep1 <- !is.na(z1) & nzchar(z1)
  sd2_df$zone[keep1] <- z1[keep1]
}

# If zone.y exists (from chemistry), fill any remaining NAs/empties
if ("zone.y" %in% names(sd2_df)) {
  z2 <- as.character(sd2_df$zone.y)   # coerce to character to avoid type clashes
  need2 <- is.na(sd2_df$zone) | !nzchar(sd2_df$zone)
  sd2_df$zone[need2] <- z2[need2]
}

# Make zone a factor for downstream models
sd2_df$zone <- factor(sd2_df$zone)

# Tidy up: drop the suffix columns if present
sd2_df <- sd2_df %>% dplyr::select(-dplyr::any_of(c("zone.x","zone.y")))

# Diagnostics CSV
readr::write_csv(
  sd2_df %>% dplyr::transmute(SampleID, site_abbrev, zone, has_chem = !is.na(contam_pc1)),
  "exports/join_diagnostics_samples_vs_chem.csv"
)

# Assign back to phyloseq
sd2_ps <- sd2_df
rownames(sd2_ps) <- sd2_ps$SampleID
sample_data(ps) <- phyloseq::sample_data(sd2_ps)

message("zone source → ",
        if ("zone.x" %in% names(sd2_df)) "metadata " else "",
        if ("zone.y" %in% names(sd2_df)) "+ chemistry" else "")









# ---------------------------
# 6) Alpha diversity
# ---------------------------
# ---------------------------
# 6) Alpha diversity
# ---------------------------
# ---------------------------
# 6) Alpha diversity  (robust SampleID handling)
# ---------------------------

# Get plain data.frame of sample metadata
sd_alpha_df <- as(phyloseq::sample_data(ps), "data.frame")

# Add a fallback rowname column (won't collide with existing SampleID)
sd_alpha_df$SampleID_rowname <- rownames(sd_alpha_df)

# Use existing SampleID if present; otherwise use rownames
if ("SampleID" %in% names(sd_alpha_df)) {
  sd_alpha_df$SampleID <- clean_id(sd_alpha_df$SampleID)
} else {
  sd_alpha_df$SampleID <- clean_id(sd_alpha_df$SampleID_rowname)
}

# Keep only the fields we need (and drop the helper column)
sd_alpha_df <- sd_alpha_df |>
  dplyr::select(
    SampleID,
    dplyr::any_of(c("site_abbrev","zone","contam_pc1","contam_index_zmean","contam_tier"))
  )

# Compute alpha and join metadata
alpha_df <- phyloseq::estimate_richness(ps, measures = c("Observed","Chao1","Shannon")) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::left_join(sd_alpha_df, by = "SampleID")

readr::write_csv(alpha_df, "exports/alpha_metrics_merged.csv")

# Rows complete for each analysis
alpha_cc_tier <- alpha_df |> dplyr::filter(!is.na(contam_tier))
alpha_cc_pc1  <- alpha_df |> dplyr::filter(!is.na(contam_pc1))

# Plots
g_alpha_box <- function(metric){
  ggplot(alpha_cc_tier, aes(x = contam_tier, y = .data[[metric]], fill = contam_tier)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_jitter(width = 0.08, height = 0, alpha = 0.7, size = 1) +
    labs(x = "Contamination tier (PC1 tertiles)", y = metric,
         title = glue::glue("{metric} by contamination tier")) +
    guides(fill="none")
}
ggsave("exports/alpha_Shannon_by_tier.png", g_alpha_box("Shannon"), width=6, height=4, dpi=300)
ggsave("exports/alpha_Chao1_by_tier.png",   g_alpha_box("Chao1"),   width=6, height=4, dpi=300)
ggsave("exports/alpha_Observed_by_tier.png",g_alpha_box("Observed"),width=6, height=4, dpi=300)

# Stats
alpha_tests <- dplyr::bind_rows(
  rstatix::kruskal_test(alpha_cc_tier, Shannon ~ contam_tier) |> dplyr::mutate(metric="Shannon"),
  rstatix::kruskal_test(alpha_cc_tier, Chao1   ~ contam_tier) |> dplyr::mutate(metric="Chao1"),
  rstatix::kruskal_test(alpha_cc_tier, Observed~ contam_tier) |> dplyr::mutate(metric="Observed")
)
readr::write_csv(alpha_tests, "exports/alpha_kruskal_by_tier.csv")

dunn_tbl <- list(
  Shannon = FSA::dunnTest(Shannon ~ contam_tier, data = alpha_cc_tier, method = "bh")$res,
  Chao1   = FSA::dunnTest(Chao1   ~ contam_tier, data = alpha_cc_tier, method = "bh")$res,
  Observed= FSA::dunnTest(Observed~ contam_tier, data = alpha_cc_tier, method = "bh")$res
) |> purrr::imap(~dplyr::mutate(.x, metric = .y)) |> dplyr::bind_rows()
readr::write_csv(dunn_tbl, "exports/alpha_dunn_by_tier.csv")

# Linear models vs PC1
# ---------------------------
# Linear models vs PC1 (robust to single-level zone)
# ---------------------------

safe_lm_alpha <- function(metric){
  d <- alpha_cc_pc1 |>
    dplyr::select(y = dplyr::all_of(metric), contam_pc1, zone) |>
    tidyr::drop_na(y, contam_pc1)
  
  # ensure factor, drop empty levels
  d$zone <- droplevels(factor(d$zone))
  has_zone <- nlevels(d$zone) >= 2
  has_var  <- stats::var(d$contam_pc1, na.rm = TRUE) > 0
  enough_n <- nrow(d) >= 3
  
  if (!enough_n || !has_var) {
    return(
      tibble::tibble(
        term      = c("(Intercept)", "contam_pc1"),
        estimate  = NA_real_, std.error = NA_real_,
        statistic = NA_real_, p.value  = NA_real_,
        metric    = metric,
        model     = if (has_zone) "y ~ contam_pc1 + zone" else "y ~ contam_pc1",
        note      = dplyr::case_when(
          !enough_n ~ "Skipped: <3 complete rows after NA drop",
          !has_var  ~ "Skipped: contam_pc1 has zero variance",
          TRUE ~ "Skipped"
        )
      )
    )
  }
  
  form <- if (has_zone) y ~ contam_pc1 + zone else y ~ contam_pc1
  fit  <- stats::lm(form, data = d)
  
  broom::tidy(fit) |>
    dplyr::mutate(
      metric = metric,
      model  = if (has_zone) "y ~ contam_pc1 + zone" else "y ~ contam_pc1",
      note   = dplyr::if_else(has_zone, NA_character_, "zone had <2 levels; dropped")
    )
}

lm_out <- dplyr::bind_rows(
  safe_lm_alpha("Shannon"),
  safe_lm_alpha("Chao1"),
  safe_lm_alpha("Observed")
)

readr::write_csv(lm_out, "exports/alpha_lm_vs_pc1_zone.csv")

# Scatter plots
g_alpha_scatter <- function(metric){
  ggplot(alpha_cc_pc1, aes(x = contam_pc1, y = .data[[metric]], color = zone)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(x = "Contamination gradient (PC1)", y = metric,
         title = glue::glue("{metric} ~ contamination (PC1)")) +
    theme(legend.position = "bottom")
}
ggsave("exports/alpha_Shannon_vs_pc1.png", g_alpha_scatter("Shannon"), width=6, height=4, dpi=300)
ggsave("exports/alpha_Chao1_vs_pc1.png",   g_alpha_scatter("Chao1"),   width=6, height=4, dpi=300)


# ---------------------------
# 7) Beta diversity: PERMANOVA, betadisper, dbRDA, envfit
# ---------------------------
# ---------------------------
# 7) Beta diversity: PERMANOVA, betadisper, dbRDA, envfit (robust)
# ---------------------------

# ---- Helpers (scoped to this section) ----
run_adonis_safe <- function(dist_mat, data_df, formula_str, outfile){
  ok_n <- nrow(data_df) >= 3 && nrow(dist_mat) >= 3
  if (!ok_n) {
    writeLines(glue::glue("PERMANOVA skipped: need >=3 samples (have {nrow(data_df)})."),
               outfile)
    return(invisible(NULL))
  }
  fit <- tryCatch(
    vegan::adonis2(as.formula(formula_str), data = data_df, permutations = 9999, by = "margin"),
    error = function(e) e
  )
  sink(outfile)
  cat("PERMANOVA (Bray)\n")
  cat("Model: ", formula_str, "\n", sep = "")
  if (inherits(fit, "error")) {
    cat("Skipped with error: ", conditionMessage(fit), "\n", sep = "")
  } else {
    print(fit)
  }
  sink()
}

run_betadisper_safe <- function(dist_mat, groups, outfile, label="group"){
  df <- data.frame(g = droplevels(groups))
  tbl <- table(df$g)
  enough <- length(tbl[tbl >= 2]) >= 2 && nrow(df) >= 4 && nrow(dist_mat) >= 4
  if (!enough) {
    msg <- paste0(
      "Betadisper skipped (", label, "). Need >=2 groups with >=2 samples.\n",
      "Group sizes: ", paste(paste(names(tbl), tbl, sep=":"), collapse=", "), "\n"
    )
    writeLines(msg, outfile)
    return(invisible(NULL))
  }
  bd <- tryCatch(betadisper(as.dist(dist_mat), df$g), error = function(e) e)
  sink(outfile)
  cat("Betadisper by ", label, "\n", sep = "")
  if (inherits(bd, "error")) {
    cat("Skipped with error: ", conditionMessage(bd), "\n", sep = "")
  } else {
    print(anova(bd))
  }
  sink()
}

# ---- Distances & ordination ----
ps_rel    <- transform_sample_counts(ps, function(x) if (sum(x)==0) x else x / sum(x))
dist_bray <- phyloseq::distance(ps_rel, method = "bray")
ord_pcoa  <- ordinate(ps_rel, method = "PCoA", distance = dist_bray)

# ---- Metadata frame & types ----
mdf <- as(sample_data(ps_rel), "data.frame")
mdf$SampleID     <- rownames(mdf)
mdf$zone         <- droplevels(factor(mdf$zone))
mdf$site_abbrev  <- droplevels(factor(mdf$site_abbrev))
mdf$contam_tier  <- droplevels(factor(mdf$contam_tier))
mdf$contam_pc1   <- suppressWarnings(as.numeric(mdf$contam_pc1))

# ---------- Overall PERMANOVA (chooses best feasible model) ----------
keep_all  <- !is.na(mdf$contam_pc1)
mdf_all   <- mdf[keep_all, , drop = FALSE]
dist_all  <- as.matrix(dist_bray)[keep_all, keep_all, drop = FALSE]

has_zone  <- nlevels(mdf_all$zone)        >= 2
has_site  <- nlevels(mdf_all$site_abbrev) >= 2
has_var   <- is.finite(var(mdf_all$contam_pc1, na.rm=TRUE)) && var(mdf_all$contam_pc1, na.rm=TRUE) > 0

form_str <- if (has_zone && has_site) {
  "as.dist(dist_all) ~ contam_pc1 + zone + site_abbrev"
} else if (has_zone) {
  "as.dist(dist_all) ~ contam_pc1 + zone"
} else if (has_site) {
  "as.dist(dist_all) ~ contam_pc1 + site_abbrev"
} else {
  "as.dist(dist_all) ~ contam_pc1"
}

if (!has_var) {
  writeLines("PERMANOVA skipped: contam_pc1 variance == 0.",
             "exports/permanova_overall.txt")
} else {
  run_adonis_safe(dist_all, mdf_all, form_str, "exports/permanova_overall.txt")
}

# ---------- Dispersion tests ----------
run_betadisper_safe(dist_all, mdf_all$zone,        "exports/betadisper_by_zone.txt",     label="zone")
run_betadisper_safe(dist_all, mdf_all$site_abbrev, "exports/betadisper_by_location.txt", label="site_abbrev")

# ---------- dbRDA (capscale) ----------
if (nrow(mdf_all) >= 3 && nrow(dist_all) >= 3 && has_var) {
  form_cap <- if (has_zone && has_site) {
    "as.dist(dist_all) ~ contam_pc1 + zone + site_abbrev"
  } else if (has_zone) {
    "as.dist(dist_all) ~ contam_pc1 + zone"
  } else if (has_site) {
    "as.dist(dist_all) ~ contam_pc1 + site_abbrev"
  } else {
    "as.dist(dist_all) ~ contam_pc1"
  }
  cap <- tryCatch(capscale(as.formula(form_cap), data = mdf_all, add = TRUE),
                  error = function(e) e)
  if (!inherits(cap, "error")) {
    saveRDS(cap, "exports/capscale_bray_model.rds")
    png("exports/ordination_dbRDA_bray.png", width=1600, height=1200, res=200)
    plot(cap, display = c("sites","bp"),
         main = paste0("dbRDA (Bray) ~ ", gsub("as.dist\\(dist_all\\) ~ ", "", form_cap)))
    dev.off()
  } else {
    writeLines(paste0("dbRDA skipped: ", conditionMessage(cap)),
               "exports/ordination_dbRDA_bray.txt")
  }
} else {
  writeLines("dbRDA skipped: insufficient usable samples or zero variance in contam_pc1.",
             "exports/ordination_dbRDA_bray.txt")
}

# ---------- envfit on unconstrained PCoA ----------
env_vars <- mdf %>%
  dplyr::select(cu, ni, zn, co, contam_pc1) %>%
  dplyr::mutate(dplyr::across(everything(), ~ suppressWarnings(as.numeric(.x))))
ef <- tryCatch(vegan::envfit(ord_pcoa, env_vars, permutations = 9999),
               error = function(e) e)
sink("exports/envfit_pcoa_bray.txt")
if (inherits(ef, "error")) {
  cat("envfit skipped: ", conditionMessage(ef), "\n", sep = "")
} else {
  print(ef)
}
sink()

# ---------- Per-location PERMANOVA (subset within each site_abbrev) ----------
dir.create("exports/per_location", showWarnings = FALSE)
sites <- levels(mdf$site_abbrev)
per_site_summary <- list()

for (s in sites) {
  idx <- which(mdf$site_abbrev == s & !is.na(mdf$contam_pc1))
  if (length(idx) < 3) {
    writeLines(glue::glue("Site {s}: skipped (n={length(idx)})."),
               file.path("exports/per_location", paste0("permanova_site_", s, ".txt")))
    next
  }
  mdf_s  <- mdf[idx, , drop = FALSE]
  dist_s <- as.matrix(dist_bray)[idx, idx, drop = FALSE]
  
  has_zone_s <- nlevels(droplevels(mdf_s$zone)) >= 2
  has_var_s  <- is.finite(var(mdf_s$contam_pc1, na.rm=TRUE)) && var(mdf_s$contam_pc1, na.rm=TRUE) > 0
  
  if (!has_var_s) {
    writeLines(glue::glue("Site {s}: PERMANOVA skipped (contam_pc1 variance == 0)."),
               file.path("exports/per_location", paste0("permanova_site_", s, ".txt")))
    next
  }
  
  form_s   <- if (has_zone_s) "as.dist(dist_s) ~ contam_pc1 + zone" else "as.dist(dist_s) ~ contam_pc1"
  out_path <- file.path("exports/per_location", paste0("permanova_site_", s, ".txt"))
  run_adonis_safe(dist_s, mdf_s, form_s, out_path)
  
  per_site_summary[[s]] <- data.frame(
    site_abbrev = s,
    n = nrow(mdf_s),
    zone_levels = nlevels(droplevels(mdf_s$zone)),
    used_model  = gsub("as.dist\\(dist_s\\) ~ ", "", form_s),
    stringsAsFactors = FALSE
  )
}
if (length(per_site_summary)) {
  per_site_summary <- dplyr::bind_rows(per_site_summary)
  readr::write_csv(per_site_summary, "exports/per_location/per_site_summary.csv")
}

# ---------- PCoA plots (by zone / by location / faceted) ----------
scores <- as.data.frame(ord_pcoa$vectors)
scores$SampleID <- rownames(scores)
scores <- scores %>%
  dplyr::left_join(mdf %>% dplyr::select(SampleID, zone, site_abbrev, contam_tier), by="SampleID")

pc1 <- round(ord_pcoa$values$Relative_eig[1]*100,1)
pc2 <- round(ord_pcoa$values$Relative_eig[2]*100,1)

# Color by zone
p_zone <- ggplot(scores, aes(Axis.1, Axis.2, color = zone)) +
  geom_point(size = 2, alpha = 0.9) +
  stat_ellipse(level = 0.68, linetype = 2, linewidth = 0.7, alpha = 0.6, na.rm = TRUE) +
  labs(x = glue::glue("PCoA1 ({pc1}%)"), y = glue::glue("PCoA2 ({pc2}%)"),
       title = "PCoA (Bray-Curtis) colored by zone") +
  theme(legend.position = "bottom")
ggsave("exports/beta_pcoa_by_zone.png", p_zone, width=7, height=5, dpi=300)

# Color by location
p_loc <- ggplot(scores, aes(Axis.1, Axis.2, color = site_abbrev)) +
  geom_point(size = 2, alpha = 0.9) +
  stat_ellipse(level = 0.68, linetype = 2, linewidth = 0.7, alpha = 0.6, na.rm = TRUE) +
  labs(x = glue::glue("PCoA1 ({pc1}%)"), y = glue::glue("PCoA2 ({pc2}%)"),
       title = "PCoA (Bray-Curtis) colored by location") +
  theme(legend.position = "bottom")
ggsave("exports/beta_pcoa_by_location.png", p_loc, width=7, height=5, dpi=300)

# Facet by zone
p_facet_zone <- ggplot(scores, aes(Axis.1, Axis.2, color = site_abbrev)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(x = glue::glue("PCoA1 ({pc1}%)"), y = glue::glue("PCoA2 ({pc2}%)"),
       title = "PCoA (Bray-Curtis) faceted by zone") +
  theme(legend.position = "bottom") +
  facet_wrap(~ zone, drop = TRUE)
ggsave("exports/beta_pcoa_facet_by_zone.png", p_facet_zone, width=9, height=6, dpi=300)

# Facet by location
p_facet_loc <- ggplot(scores, aes(Axis.1, Axis.2, color = zone)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(x = glue::glue("PCoA1 ({pc1}%)"), y = glue::glue("PCoA2 ({pc2}%)"),
       title = "PCoA (Bray-Curtis) faceted by location") +
  theme(legend.position = "bottom") +
  facet_wrap(~ site_abbrev, drop = TRUE)
ggsave("exports/beta_pcoa_facet_by_location.png", p_facet_loc, width=10, height=7, dpi=300)

# Export scores
readr::write_csv(scores, "exports/pcoa_scores_with_zone_and_location.csv")


# ---------------------------
# 8) Taxa ~ contaminants (correlations + optional ANCOM-BC2)
# ---------------------------
# ---------------------------
# 8) Taxa ~ contaminants (correlations + optional ANCOM-BC2) — ROBUST
# ---------------------------

# Helper: pick deepest rank with enough non-NA/non-empty labels
pick_rank <- function(ps, candidates = c("Genus","genus","Family","family","Order","order","Class","class","Phylum","phylum","Kingdom","kingdom")) {
  if (is.null(tax_table(ps, errorIfNULL = FALSE))) return(NA_character_)
  tt <- as.data.frame(tax_table(ps))
  for (rk in candidates) {
    if (rk %in% colnames(tt)) {
      vals <- as.character(tt[[rk]])
      nn <- sum(!is.na(vals) & nzchar(vals) & vals != "NA")
      if (nn >= 2) return(rk)
    }
  }
  NA_character_
}

collapse_safe <- function(ps, rk) {
  if (is.na(rk)) return(list(obj = ps, used = "ASV (no collapse; no usable rank)"))
  # keep NA-labeled taxa together instead of dropping them
  tmp <- try(suppressWarnings(phyloseq::tax_glom(ps, taxrank = rk, NArm = FALSE)), silent = TRUE)
  if (inherits(tmp, "try-error") || ntaxa(tmp) == 0) {
    list(obj = ps, used = paste0("ASV (fallback; collapse at ", rk, " failed)"))
  } else {
    list(obj = tmp, used = rk)
  }
}

rk <- pick_rank(ps)
col_res <- collapse_safe(ps, rk)
ps_coll <- col_res$obj
used_rank <- col_res$used
message("Taxa collapse level used: ", used_rank)

# Relative abundance for correlations
ps_coll_rel <- transform_sample_counts(ps_coll, function(x) if (sum(x)==0) x else x / sum(x))

# Pull tables
otu_df <- as.data.frame(otu_table(ps_coll_rel))
tax_df <- if (!is.null(tax_table(ps_coll_rel, errorIfNULL = FALSE))) {
  as.data.frame(tax_table(ps_coll_rel)) %>% tibble::rownames_to_column("Feature")
} else {
  tibble::tibble(Feature = rownames(otu_df))
}

# Long table: SampleID x Feature x rel_abund
long_df <- otu_df %>%
  tibble::rownames_to_column("Feature") %>%
  tidyr::pivot_longer(-Feature, names_to = "SampleID", values_to = "rel_abund")

# Choose a display name column for taxa
tax_disp_col <- intersect(c("Genus","genus","Family","family","Order","order","Class","class","Phylum","phylum","Kingdom","kingdom"), colnames(tax_df))
if (length(tax_disp_col) == 0) tax_disp_col <- "Feature"
disp_name <- if (length(tax_disp_col) > 0) tax_disp_col[1] else "Feature"

long_df <- long_df %>%
  dplyr::left_join(tax_df %>% dplyr::select(Feature, !!disp_name), by = "Feature") %>%
  dplyr::rename(Taxon = !!disp_name)

# Numerical env variables from the Section 7 metadata (mdf)
stopifnot(exists("mdf"))
num_vars <- c("cu","ni","zn","co","contam_pc1")
mdf_num <- mdf %>%
  dplyr::transmute(SampleID = rownames(.),
                   cu = suppressWarnings(as.numeric(cu)),
                   ni = suppressWarnings(as.numeric(ni)),
                   zn = suppressWarnings(as.numeric(zn)),
                   co = suppressWarnings(as.numeric(co)),
                   contam_pc1 = suppressWarnings(as.numeric(contam_pc1)))

# Spearman correlations Taxon vs metals/PC1 (per taxon), BH FDR
safe_cor_test <- function(x, y) {
  # Handle constant vectors or all-NA gracefully
  if (all(is.na(x)) || all(is.na(y))) return(c(rho = NA_real_, p = NA_real_))
  if (stats::var(x, na.rm = TRUE) == 0 || stats::var(y, na.rm = TRUE) == 0) return(c(rho = NA_real_, p = NA_real_))
  ct <- try(stats::cor.test(x, y, method = "spearman"), silent = TRUE)
  if (inherits(ct, "try-error")) return(c(rho = NA_real_, p = NA_real_))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

gen_cor <- long_df %>%
  dplyr::filter(!is.na(Taxon) & nzchar(Taxon)) %>%
  dplyr::left_join(mdf_num, by = "SampleID") %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarize(
    # rho
    rho_cu = safe_cor_test(rel_abund, cu)["rho"],
    rho_ni = safe_cor_test(rel_abund, ni)["rho"],
    rho_zn = safe_cor_test(rel_abund, zn)["rho"],
    rho_co = safe_cor_test(rel_abund, co)["rho"],
    rho_contam_pc1 = safe_cor_test(rel_abund, contam_pc1)["rho"],
    # p
    p_cu = safe_cor_test(rel_abund, cu)["p"],
    p_ni = safe_cor_test(rel_abund, ni)["p"],
    p_zn = safe_cor_test(rel_abund, zn)["p"],
    p_co = safe_cor_test(rel_abund, co)["p"],
    p_contam_pc1 = safe_cor_test(rel_abund, contam_pc1)["p"],
    .groups = "drop"
  ) %>%
  # FDR for each p-column
  dplyr::mutate(
    q_cu = p.adjust(p_cu, method = "BH"),
    q_ni = p.adjust(p_ni, method = "BH"),
    q_zn = p.adjust(p_zn, method = "BH"),
    q_co = p.adjust(p_co, method = "BH"),
    q_contam_pc1 = p.adjust(p_contam_pc1, method = "BH"),
    collapse_level = used_rank
  ) %>%
  # For safety, coerce to numeric
  dplyr::mutate(dplyr::across(dplyr::starts_with(c("rho_","p_","q_")), ~ as.numeric(.)))

readr::write_csv(gen_cor, "exports/genus_spearman_correlations_vs_metals_pc1.csv")

# Heatmap of top taxa vs contaminants (if enough)
sig_taxa <- gen_cor %>%
  dplyr::arrange(q_contam_pc1) %>%
  dplyr::filter(!is.na(q_contam_pc1)) %>%
  dplyr::slice_head(n = 40)

if (nrow(sig_taxa) >= 2) {
  mat <- sig_taxa %>%
    dplyr::select(rho_cu, rho_ni, rho_zn, rho_co, rho_contam_pc1) %>%
    as.matrix()
  rownames(mat) <- sig_taxa$Taxon
  png("exports/heatmap_taxa_rho_vs_metals_pc1.png", width=1400, height=1600, res=200)
  pheatmap::pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                     main = paste0("Taxa (", used_rank, ") ~ contaminants (Spearman rho)"))
  dev.off()
}

# ---------------------------
# Optional: ANCOM-BC2 (if installed & enough replication)
# ---------------------------
if (requireNamespace("ANCOMBC", quietly = TRUE)) {
  message("ANCOM-BC2: checking replication by contamination tier...")
  an_md <- as(phyloseq::sample_data(ps_coll), "data.frame")
  # Use contamination tier; require >=3 per group
  if ("contam_tier" %in% names(an_md)) {
    an_md$contam_tier <- factor(an_md$contam_tier, levels = c("Low","Medium","High"))
  } else {
    an_md$contam_tier <- NA
  }
  ok_rep <- all(table(an_md$contam_tier) >= 3, na.rm = TRUE) && nlevels(droplevels(an_md$contam_tier)) >= 2
  if (ok_rep) {
    out <- try(ANCOMBC::ancombc2(
      data = ps_coll,                # use collapsed counts (not relative)
      fix_formula = "contam_tier + zone",
      p_adj_method = "BH",
      prv_cut = 0.05, lib_cut = 1000,
      s0_perc = 0.05, alpha = 0.1, n_cl = 1
    ), silent = TRUE)
    if (!inherits(out, "try-error")) {
      tx_tbl <- as.data.frame(tax_table(ps_coll)) %>%
        tibble::rownames_to_column("Feature")
      da_tab <- out$res$diff_abn %>%
        tibble::rownames_to_column("Feature") %>%
        dplyr::left_join(tx_tbl, by = "Feature") %>%
        dplyr::arrange(q_val)
      readr::write_csv(da_tab, "exports/ancombc2_taxa_contam_tier.csv")
    } else {
      writeLines(paste0("ANCOM-BC2 skipped: ", as.character(out)),
                 "exports/ancombc2_taxa_contam_tier.txt")
    }
  } else {
    writeLines("ANCOM-BC2 skipped: need ≥3 samples per contamination tier (and ≥2 tiers present).",
               "exports/ancombc2_taxa_contam_tier.txt")
  }
} else {
  writeLines("ANCOM-BC2 not installed; skipping differential abundance.",
             "exports/ancombc2_taxa_contam_tier.txt")
}


# ---------------------------
# 9) Session info & README
# ---------------------------
readr::write_lines(c(
  "Vieques Microbiome Pipeline - Outputs overview",
  "exports/alpha_* : alpha diversity plots & stats",
  "exports/permanova_* : PERMANOVA results",
  "exports/betadisper_* : dispersion tests",
  "exports/ordination_* : ordination figures",
  "exports/capscale_* : dbRDA object",
  "exports/envfit_* : envfit loadings",
  "exports/genus_* : genus-level correlation/DA results",
  "exports/pcoa_scores_with_groups.csv : ordination scores table",
  "exports/join_diagnostics_samples_vs_chem.csv : sample↔chem join success table",
  "exports/chem_clean_with_index.csv : chemistry + contamination indices used"
), "exports/README_outputs.txt")

sink("exports/sessionInfo.txt"); print(sessionInfo()); sink()
message("✅ Pipeline complete. See the 'exports/' folder for results.")





# ===============================================================
# view_exports.R — quick viewer + HTML gallery for pipeline outputs
# ===============================================================

# Packages (only common, light deps)
need <- c("readr","dplyr","tibble","stringr","janitor","png","grid")
for(p in need){
  if(!requireNamespace(p, quietly = TRUE)){
    install.packages(p, repos="https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

exp_dir <- "exports"
stopifnot(dir.exists(exp_dir))

# Gather files
all_files <- list.files(exp_dir, full.names = TRUE)
pngs <- all_files[str_detect(tolower(all_files), "\\.png$")]
csvs <- all_files[str_detect(tolower(all_files), "\\.csv$")]
txts <- all_files[str_detect(tolower(all_files), "\\.txt$")]
rds  <- all_files[str_detect(tolower(all_files), "\\.rds$")]
others <- setdiff(all_files, c(pngs, csvs, txts, rds))

# ---------------------------
# 1) Console summary
# ---------------------------
cat("\n=== EXPORTS SUMMARY ===\n")
cat("PNG plots: ", length(pngs), "\n")
cat("CSV tables:", length(csvs), "\n")
cat("TXT files: ", length(txts), "\n")
cat("RDS files: ", length(rds), "\n")
cat("Other:     ", length(others), "\n\n")

# ---------------------------
# 2) Optional: show all PNGs in R plotting window
#     (comment this block if you don't want pop-ups)
# ---------------------------
if(length(pngs)){
  cat("Showing PNGs in the plotting window...\n")
  for (p in sort(pngs)) {
    try({
      img <- png::readPNG(p)
      grid::grid.newpage(); grid::grid.raster(img)
      title_txt <- basename(p)
      # draw a lightweight title
      grid::grid.text(title_txt, x=0.02, y=0.98, just=c("left","top"))
      Sys.sleep(0.5)
    }, silent = TRUE)
  }
}

# ---------------------------
# 3) Build an HTML gallery (exports/_index.html)
# ---------------------------
# helper: html escape
html_escape <- function(x){
  x <- gsub("&","&amp;",x, fixed=TRUE)
  x <- gsub("<","&lt;",x, fixed=TRUE)
  x <- gsub(">","&gt;",x, fixed=TRUE)
  x
}

# read first n lines of a text file
head_text <- function(path, n = 40){
  if(!file.exists(path)) return("")
  ln <- try(readLines(path, warn = FALSE), silent = TRUE)
  if(inherits(ln, "try-error")) return("")
  paste(utils::head(ln, n), collapse = "\n")
}

# read first rows of a CSV safely
head_csv_html <- function(path, n = 10){
  out <- try(suppressMessages(readr::read_csv(path, n_max = n, show_col_types = FALSE)), silent = TRUE)
  if(inherits(out, "try-error")) return("<em>(Could not preview CSV)</em>")
  out <- janitor::remove_empty(out, c("rows","cols"))
  # build simple HTML table
  if(nrow(out)==0) return("<em>(Empty CSV)</em>")
  hdr <- paste(sprintf("<th>%s</th>", html_escape(colnames(out))), collapse = "")
  rows <- apply(out, 1, function(r){
    cells <- paste(sprintf("<td>%s</td>", html_escape(as.character(r))), collapse = "")
    sprintf("<tr>%s</tr>", cells)
  })
  paste0("<table class='mini'><thead><tr>", hdr, "</tr></thead><tbody>",
         paste(rows, collapse = "\n"), "</tbody></table>")
}

# Minimal CSS
css <- "
body { font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif; margin: 24px; }
h1 { margin-top: 0; }
.grid { display: grid; gap: 16px; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); }
.card { border: 1px solid #e5e7eb; border-radius: 12px; padding: 12px; box-shadow: 0 1px 2px rgba(0,0,0,.03); }
.card img { max-width: 100%; height: auto; border-radius: 8px; display: block; }
.kv { color: #555; font-size: 0.9em; margin: 0 0 6px 0; }
a { color: #0b70ff; text-decoration: none; }
a:hover { text-decoration: underline; }
.mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; white-space: pre-wrap; background: #f8fafc; padding: 8px; border-radius: 8px; border: 1px solid #e5e7eb; }
.mini { border-collapse: collapse; width: 100%; font-size: 12px; }
.mini th, .mini td { border: 1px solid #e5e7eb; padding: 4px 6px; vertical-align: top; }
.section { margin-top: 28px; }
footer { margin-top: 24px; color: #666; font-size: 12px; }
"

html <- c(
  "<!doctype html>",
  "<html><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width, initial-scale=1'>",
  "<title>Exports Gallery</title>",
  "<style>", css, "</style>",
  "</head><body>",
  "<h1>Exports Gallery</h1>",
  sprintf("<p class='kv'>Directory: <span class='mono'>%s</span></p>", html_escape(normalizePath(exp_dir)))
)

# PNG section
html <- c(html, "<div class='section'><h2>Plots (PNG)</h2>")
if(length(pngs)){
  html <- c(html, "<div class='grid'>")
  for (p in sort(pngs)) {
    rel <- file.path(basename(exp_dir), basename(p))
    html <- c(html,
              "<div class='card'>",
              sprintf("<div class='kv'>%s</div>", html_escape(basename(p))),
              sprintf("<a href='%s'><img src='%s' alt='%s'></a>", html_escape(rel), html_escape(rel), html_escape(basename(p))),
              "</div>")
  }
  html <- c(html, "</div>")
} else {
  html <- c(html, "<p><em>No PNG plots found.</em></p>")
}
html <- c(html, "</div>")

# CSV section
html <- c(html, "<div class='section'><h2>Tables (CSV)</h2>")
if(length(csvs)){
  for (cpath in sort(csvs)) {
    rel <- file.path(basename(exp_dir), basename(cpath))
    html <- c(html,
              "<div class='card'>",
              sprintf("<div class='kv'><a href='%s'>%s</a></div>", html_escape(rel), html_escape(basename(cpath))),
              head_csv_html(cpath),
              "</div>")
  }
} else {
  html <- c(html, "<p><em>No CSV files found.</em></p>")
}
html <- c(html, "</div>")

# TXT section
html <- c(html, "<div class='section'><h2>Text Results (TXT)</h2>")
if(length(txts)){
  for (tpath in sort(txts)) {
    rel <- file.path(basename(exp_dir), basename(tpath))
    html <- c(html,
              "<div class='card'>",
              sprintf("<div class='kv'><a href='%s'>%s</a></div>", html_escape(rel), html_escape(basename(tpath))),
              sprintf("<div class='mono'>%s</div>", html_escape(head_text(tpath, 80))),
              "</div>")
  }
} else {
  html <- c(html, "<p><em>No TXT files found.</em></p>")
}
html <- c(html, "</div>")

# RDS section
html <- c(html, "<div class='section'><h2>Objects (RDS)</h2>")
if(length(rds)){
  html <- c(html, "<ul>")
  for (r in sort(rds)) {
    rel <- file.path(basename(exp_dir), basename(r))
    html <- c(html, sprintf("<li><a href='%s'>%s</a> — <span class='kv'>R object (not previewed)</span></li>",
                            html_escape(rel), html_escape(basename(r))))
  }
  html <- c(html, "</ul>")
} else {
  html <- c(html, "<p><em>No RDS files found.</em></p>")
}
html <- c(html, "</div>")

# Other files
if(length(others)){
  html <- c(html, "<div class='section'><h2>Other files</h2><ul>")
  for (o in sort(others)) {
    rel <- file.path(basename(exp_dir), basename(o))
    html <- c(html, sprintf("<li><a href='%s'>%s</a></li>", html_escape(rel), html_escape(basename(o))))
  }
  html <- c(html, "</ul></div>")
}

html <- c(html,
          "<footer>Generated by view_exports.R</footer>",
          "</body></html>"
)

# Write HTML
out_html <- file.path(exp_dir, "_index.html")
writeLines(html, out_html)
cat("✅ Wrote HTML gallery:", out_html, "\n")

# Open in default browser (comment out if undesired)
try(utils::browseURL(out_html), silent = TRUE)







# ===============================================================
# PCA of contaminants + Genus centroids overlay
# Produces: 
#  - exports/pca_contaminants_samples_w_genus_centroids.png
#  - exports/genus_centroids_on_contaminants_pca.csv
#  - exports/pca_contaminants_scores_by_sample.csv
#  - exports/pca_contaminants_loadings.csv
# ===============================================================

# ===============================================================
# PCA of contaminants + organics with OVERALL & BY-ZONE Genus centroids
# Outputs (saved in exports/):
#  - chemorg_pca_scores_by_site.csv
#  - chemorg_pca_loadings.csv
#  - chemorg_pca_samples.png
#  - chemorg_pca_samples_w_genus_centroids.png
#  - chemorg_pca_genus_centroids_overall.csv
#  - chemorg_pca_genus_centroids_by_zone.csv
#  - chemorg_pca_genus_centroids_by_zone.png
#  - chemorg_pca_biplot_loadings.png
# ===============================================================

# ---- deps# ===============================================================
# Flexible PCA plotting: by ZONE, by LOCATION, or SAMPLE ID
# Also: genus centroids overall / by zone / by location
# Files written under exports/
# ===============================================================

# --- helpers ---
.ensure_pkg <- function(p){ if(!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org"); suppressPackageStartupMessages(library(p, character.only=TRUE)) }
.ensure_pkg("ggrepel")

# ensures factor with >=2 levels to use as color; else returns NULL
.valid_group <- function(df, col){
  if (!col %in% names(df)) return(NULL)
  f <- droplevels(as.factor(df[[col]]))
  if (nlevels(f) >= 2) f else NULL
}

# Build a sample plot in different modes: "zone", "location", "samples"
plot_chem_pca <- function(sample_scores, pc1v, pc2v, mode = c("zone","location","samples"), outfile){
  mode <- match.arg(mode)
  df <- sample_scores
  has_zone <- !is.null(.valid_group(df, "zone"))
  has_loc  <- !is.null(.valid_group(df, "site_abbrev"))
  
  p <- ggplot(df, aes(PC1, PC2, label = SampleID))
  
  if (mode == "zone" && has_zone){
    p <- p + aes(color = zone, shape = site_abbrev)
    leg_sub <- c(color = "Zone", shape = "Location")
  } else if (mode == "location" && has_loc){
    p <- p + aes(color = site_abbrev)
    leg_sub <- c(color = "Location")
  } else {
    leg_sub <- NULL
  }
  
  p <- p +
    geom_point(size = 2.2, alpha = 0.9, na.rm = TRUE) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 60, show.legend = FALSE) +
    labs(
      title = paste("Contaminants + Organics PCA —", switch(mode, zone="by Zone", location="by Location", samples="per Sample")),
      x = glue::glue("PC1 ({pc1v}%)"), y = glue::glue("PC2 ({pc2v}%)")
    ) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  
  if (!is.null(leg_sub)){
    if ("color" %in% names(leg_sub) && "shape" %in% names(leg_sub)){
      p <- p + labs(color = leg_sub[["color"]], shape = leg_sub[["shape"]])
    } else if ("color" %in% names(leg_sub)){
      p <- p + labs(color = leg_sub[["color"]])
    }
  }
  
  ggsave(outfile, p, width = 9, height = 7, dpi = 300)
  p
}

# ------- SAMPLE-LEVEL PLOTS (three variants) -------
# Uses your existing sample_scores (from earlier block). 
# These will auto-fallback if zone/location is not informative.
plot_chem_pca(sample_scores, pc1v, pc2v, mode = "zone",     outfile = "exports/chemorg_pca_samples_by_zone.png")
plot_chem_pca(sample_scores, pc1v, pc2v, mode = "location", outfile = "exports/chemorg_pca_samples_by_location.png")
plot_chem_pca(sample_scores, pc1v, pc2v, mode = "samples",  outfile = "exports/chemorg_pca_samples_by_sampleid.png")

# ------- GENUS CENTROIDS: overall, by zone, by location -------
# long_abund built earlier: Feature, SampleID, rel_abund, Genus, site_abbrev, zone, PC1, PC2

# Overall centroid per Genus (weighted by rel_abund)
gen_centroids_overall <- long_abund |>
  dplyr::group_by(Genus) |>
  dplyr::summarise(
    w   = sum(rel_abund, na.rm = TRUE),
    PC1 = ifelse(w > 0, sum(PC1 * rel_abund, na.rm = TRUE) / w, NA_real_),
    PC2 = ifelse(w > 0, sum(PC2 * rel_abund, na.rm = TRUE) / w, NA_real_)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(is.finite(PC1), is.finite(PC2)) |>
  dplyr::arrange(dplyr::desc(w))
readr::write_csv(gen_centroids_overall, "exports/chemorg_pca_genus_centroids_overall.csv")

# By ZONE
gen_centroids_zone <- long_abund |>
  dplyr::group_by(Genus, zone) |>
  dplyr::summarise(
    w   = sum(rel_abund, na.rm = TRUE),
    PC1 = ifelse(w > 0, sum(PC1 * rel_abund, na.rm = TRUE) / w, NA_real_),
    PC2 = ifelse(w > 0, sum(PC2 * rel_abund, na.rm = TRUE) / w, NA_real_)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(is.finite(PC1), is.finite(PC2))
readr::write_csv(gen_centroids_zone, "exports/chemorg_pca_genus_centroids_by_zone.csv")

# By LOCATION
gen_centroids_loc <- long_abund |>
  dplyr::group_by(Genus, site_abbrev) |>
  dplyr::summarise(
    w   = sum(rel_abund, na.rm = TRUE),
    PC1 = ifelse(w > 0, sum(PC1 * rel_abund, na.rm = TRUE) / w, NA_real_),
    PC2 = ifelse(w > 0, sum(PC2 * rel_abund, na.rm = TRUE) / w, NA_real_)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(is.finite(PC1), is.finite(PC2))
readr::write_csv(gen_centroids_loc, "exports/chemorg_pca_genus_centroids_by_location.csv")




# ------- PLOTS with centroids (overall / by zone / by location) -------
topN <- 25
lab_overall <- gen_centroids_overall |> dplyr::slice_head(n = min(topN, n()))

# Overall
p_overall <- ggplot() +
  geom_point(data = sample_scores, aes(PC1, PC2, color = .valid_group(sample_scores, "zone") %||% NULL, shape = .valid_group(sample_scores, "site_abbrev") %||% NULL),
             size = 2, alpha = 0.85, na.rm = TRUE, show.legend = TRUE) +
  geom_point(data = gen_centroids_overall, aes(PC1, PC2, size = w), color = "black", alpha = 0.35, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = lab_overall, aes(PC1, PC2, label = Genus), size = 3, max.overlaps = 60, seed = 1) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  scale_size_continuous(name = "Overall genus abundance", range = c(1.2, 8)) +
  labs(title = glue::glue("Chem+Org PCA with Genus Centroids (top {nrow(lab_overall)} labeled)"),
       x = glue::glue("PC1 ({pc1v}%)"), y = glue::glue("PC2 ({pc2v}%)")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave("exports/chemorg_pca_samples_w_genus_centroids.png", p_overall, width = 9, height = 7, dpi = 300)

# By ZONE (facet)
keep_gen <- lab_overall$Genus
genz_plot <- gen_centroids_zone |> dplyr::filter(Genus %in% keep_gen)
p_by_zone <- ggplot() +
  geom_point(data = sample_scores, aes(PC1, PC2, color = zone), size = 1.8, alpha = 0.7, na.rm = TRUE) +
  geom_point(data = genz_plot, aes(PC1, PC2, size = w), color = "black", alpha = 0.45, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = genz_plot, aes(PC1, PC2, label = Genus), size = 2.8, seed = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  scale_size_continuous(name = "Genus abundance in zone", range = c(1.0, 7)) +
  labs(title = glue::glue("Genus centroids by Zone in Chem+Org PCA (top {length(keep_gen)} genera)"),
       x = glue::glue("PC1 ({pc1v}%)"), y = glue::glue("PC2 ({pc2v}%)")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  facet_wrap(~ zone, drop = TRUE)
ggsave("exports/chemorg_pca_genus_centroids_by_zone.png", p_by_zone, width = 11, height = 7.5, dpi = 300)

# By LOCATION (facet)
genl_plot <- gen_centroids_loc |> dplyr::filter(Genus %in% keep_gen)
p_by_loc <- ggplot() +
  geom_point(data = sample_scores, aes(PC1, PC2, color = site_abbrev), size = 1.8, alpha = 0.7, na.rm = TRUE) +
  geom_point(data = genl_plot, aes(PC1, PC2, size = w), color = "black", alpha = 0.45, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = genl_plot, aes(PC1, PC2, label = Genus), size = 2.5, seed = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.3, color = "grey50") +
  scale_size_continuous(name = "Genus abundance at location", range = c(1.0, 7)) +
  labs(title = glue::glue("Genus centroids by Location in Chem+Org PCA (top {length(keep_gen)} genera)"),
       x = glue::glue("PC1 ({pc1v}%)"), y = glue::glue("PC2 ({pc2v}%)")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  facet_wrap(~ site_abbrev, drop = TRUE)
ggsave("exports/chemorg_pca_genus_centroids_by_location.png", p_by_loc, width = 12, height = 8.5, dpi = 300)


message("✅ Chem+Org PCA with genus centroids complete. See files in exports/.")















######
######
######


# ===============================================================
# Phylogenetic tree + per-site heatmaps (top 30 taxa)
# ===============================================================

dir.create("exports", showWarnings = FALSE)
dir.create("exports/tree", showWarnings = FALSE)
dir.create("exports/heatmaps_by_site", showWarnings = FALSE)

# ---------- 1) Build/attach a phylogenetic tree ----------
# Use existing tree if present; otherwise create a taxonomy-based proxy tree.
if (is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))) {
  message("No phylogenetic tree detected; building a taxonomy-proxy tree.")
  tt <- as.data.frame(phyloseq::tax_table(ps))
  if (nrow(tt) < 2) stop("Need >=2 taxa to build a tree.")
  tt[is.na(tt)] <- ""
  ranks <- intersect(c("phylum","class","order","family","genus","species"), tolower(colnames(tt)))
  if (length(ranks) == 0) ranks <- colnames(tt)
  
  # distance = number of rank mismatches (simple, robust taxonomy proxy)
  n <- nrow(tt)
  D <- matrix(0, n, n, dimnames = list(rownames(tt), rownames(tt)))
  for (i in seq_len(n)) {
    for (j in i:n) {
      d <- sum(as.character(tt[i, ranks, drop=FALSE]) != as.character(tt[j, ranks, drop=FALSE]))
      D[i, j] <- D[j, i] <- d
    }
  }
  hc <- stats::hclust(as.dist(D), method = "average")
  tree_proxy <- ape::as.phylo(hc)
  phyloseq::phy_tree(ps) <- tree_proxy
}

# Save a quick tree plot of ALL taxa
png("exports/tree/tree_all_taxa.png", width = 1600, height = 1200, res = 200)
plot(phyloseq::phy_tree(ps), main = "Phylogenetic (proxy) tree — all taxa")
dev.off()

# ---------- 2) Prep: collapse to a readable taxonomic rank ----------
# Prefer Genus → Family → Order → Phylum → (ASV fallback)
tax_tbl <- as.data.frame(phyloseq::tax_table(ps))
rank_order <- c("Genus","genus","Family","family","Order","order","Phylum","phylum")
use_rank <- intersect(rank_order, colnames(tax_tbl))
use_rank <- if (length(use_rank)) use_rank[1] else NA_character_

ps_tx <- if (!is.na(use_rank)) suppressWarnings(phyloseq::tax_glom(ps, taxrank = use_rank, NArm = FALSE)) else ps
ps_tx_rel <- phyloseq::transform_sample_counts(ps_tx, function(x) if (sum(x) == 0) x else x / sum(x))

# Helper to label taxa nicely
label_for_taxa <- function(ps_obj, row_ids){
  tx <- as.data.frame(phyloseq::tax_table(ps_obj))
  if (!is.na(use_rank) && use_rank %in% names(tx)) {
    out <- as.character(tx[row_ids, use_rank])
    out[is.na(out) | !nzchar(out)] <- row_ids[is.na(out) | !nzchar(out)]
  } else {
    out <- row_ids
  }
  make.names(out, unique = TRUE)
}

# ---------- 3) Per-site heatmaps: top 30 taxa within each site ----------
sd <- as(phyloseq::sample_data(ps_tx_rel), "data.frame")
stopifnot("site_abbrev" %in% names(sd))
sites <- sort(unique(as.character(sd$site_abbrev)))

for (s in sites) {
  # subset to site
  ps_s <- phyloseq::subset_samples(ps_tx_rel, site_abbrev == s)
  ps_s <- phyloseq::prune_samples(phyloseq::sample_sums(ps_s) > 0, ps_s)
  if (phyloseq::nsamples(ps_s) < 1 || phyloseq::ntaxa(ps_s) < 2) {
    message("Skipping site ", s, " (not enough data).")
    next
  }
  
  # matrix (taxa x samples), pick top 30 by total rel. abundance IN THIS SITE
  M <- as(phyloseq::otu_table(ps_s), "matrix")
  if (!phyloseq::taxa_are_rows(ps_s)) M <- t(M)
  keep <- order(rowSums(M, na.rm = TRUE), decreasing = TRUE)[seq_len(min(30, nrow(M)))]
  M <- M[keep, , drop = FALSE]
  
  # row labels
  rownames(M) <- label_for_taxa(ps_s, rownames(M))
  
  # sample annotations (optional if these exist)
  ann <- as(phyloseq::sample_data(ps_s), "data.frame")
  ann$SampleID <- rownames(ann)
  ann_df <- ann %>%
    dplyr::transmute(
      SampleID,
      Zone = if ("zone" %in% names(ann)) ann$zone else NA,
      Organic = if ("organic_pct" %in% names(ann)) suppressWarnings(as.numeric(ann$organic_pct)) else NA_real_,
      Survival = if ("survival_pct" %in% names(ann)) suppressWarnings(as.numeric(ann$survival_pct)) else NA_real_
    ) %>% tibble::column_to_rownames("SampleID")
  
  # ensure column order aligned to annotation
  M <- M[, rownames(ann_df), drop = FALSE]
  
  # heatmap
  out_png <- file.path("exports/heatmaps_by_site", paste0("heatmap_top30_", s, ".png"))
  pheatmap::pheatmap(
    M,
    scale = "row",                        # emphasize within-taxon variation across samples
    clustering_distance_cols = "euclidean",
    clustering_method = "average",
    annotation_col = ann_df,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    main = paste0("Top 30 taxa — site ", s),
    filename = out_png,
    width = 11, height = 8
  )
}

# ---------- 4) (Optional) Tree of top 30 taxa per site ----------
# Saves a small tree per site for the same taxa shown in each heatmap.
for (s in sites) {
  ps_s <- phyloseq::subset_samples(ps_tx_rel, site_abbrev == s)
  if (phyloseq::nsamples(ps_s) < 1 || phyloseq::ntaxa(ps_s) < 2) next
  M <- as(phyloseq::otu_table(ps_s), "matrix"); if (!phyloseq::taxa_are_rows(ps_s)) M <- t(M)
  keep <- order(rowSums(M, na.rm = TRUE), decreasing = TRUE)[seq_len(min(30, nrow(M)))]
  taxa_keep <- rownames(M)[keep]
  
  # prune phyloseq to the top taxa for this site
  ps_top <- phyloseq::prune_taxa(taxa_keep, ps_s)
  
  # plot tree (proxy or real, depending on ps)
  out_tree <- file.path("exports/tree", paste0("tree_top30_", s, ".png"))
  png(out_tree, width = 1400, height = 1100, res = 200)
  plot(phyloseq::phy_tree(ps_top),
       main = paste0("Phylogenetic (proxy) tree — top 30 taxa at ", s))
  dev.off()
}

message("✅ Wrote: exports/tree/*.png and exports/heatmaps_by_site/heatmap_top30_*.png")



######
######
######
######
# ===============================================================
# Finish: tree plot + per-site heatmaps (top 30 taxa)
# (safe even if you've run parts already)
# ===============================================================

dir.create("exports", showWarnings = FALSE)
dir.create("exports/tree", showWarnings = FALSE)
dir.create("exports/heatmaps_by_site", showWarnings = FALSE)

# ---------- A) Make/attach a phylogenetic tree if missing ----------
if (is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))) {
  message("No phylogenetic tree detected; building a taxonomy-proxy tree.")
  tt <- as.data.frame(phyloseq::tax_table(ps))
  stopifnot(nrow(tt) >= 2)
  tt[is.na(tt)] <- ""
  # choose available ranks, case-insensitive
  ranks <- tolower(colnames(tt))
  ranks <- intersect(c("phylum","class","order","family","genus","species"), ranks)
  if (!length(ranks)) ranks <- tolower(colnames(tt))
  n <- nrow(tt)
  D <- matrix(0, n, n, dimnames = list(rownames(tt), rownames(tt)))
  for (i in seq_len(n)) {
    for (j in i:n) {
      d <- sum(as.character(tt[i, ranks, drop=FALSE]) != as.character(tt[j, ranks, drop=FALSE]))
      D[i, j] <- D[j, i] <- d
    }
  }
  hc <- stats::hclust(as.dist(D), method = "average")
  phyloseq::phy_tree(ps) <- ape::as.phylo(hc)
}

# Save a quick tree plot of ALL taxa (PNG)
png("exports/tree/tree_all_taxa.png", width = 1600, height = 1200, res = 200)
plot(phyloseq::phy_tree(ps), main = "Phylogenetic (proxy) tree — all taxa")
dev.off()

# ---------- B) Collapse to a readable rank & relative abundance ----------
tax_tbl <- as.data.frame(phyloseq::tax_table(ps))
rank_order <- c("Genus","genus","Family","family","Order","order","Phylum","phylum")
use_rank <- intersect(rank_order, colnames(tax_tbl))
use_rank <- if (length(use_rank)) use_rank[1] else NA_character_

ps_tx <- if (!is.na(use_rank)) suppressWarnings(phyloseq::tax_glom(ps, taxrank = use_rank, NArm = FALSE)) else ps
ps_tx_rel <- phyloseq::transform_sample_counts(ps_tx, function(x) if (sum(x) == 0) x else x / sum(x))

label_for_taxa <- function(ps_obj, row_ids){
  tx <- as.data.frame(phyloseq::tax_table(ps_obj))
  if (!is.na(use_rank) && use_rank %in% names(tx)) {
    out <- as.character(tx[row_ids, use_rank])
    out[is.na(out) | !nzchar(out)] <- row_ids[is.na(out) | !nzchar(out)]
  } else out <- row_ids
  make.names(out, unique = TRUE)
}

# ---------- C) Ensure site_abbrev exists in sample_data ----------
sd <- as(phyloseq::sample_data(ps_tx_rel), "data.frame")
if (!"site_abbrev" %in% names(sd)) {
  # try to infer from common alternatives
  guess <- intersect(names(sd), c("site","site_name","site_abbreviation","site_abbrev","site_abbr","abbrev","abbr"))
  stopifnot(length(guess) >= 1)
  sd$site_abbrev <- toupper(trimws(as.character(sd[[guess[1]]])))
  sample_data(ps_tx_rel)$site_abbrev <- sd$site_abbrev
}
sites <- sort(unique(as.character(sample_data(ps_tx_rel)$site_abbrev)))

# ---------- D) Per-site heatmaps: top 30 taxa ----------
for (s in sites) {
  ps_s <- phyloseq::subset_samples(ps_tx_rel, site_abbrev == s)
  ps_s <- phyloseq::prune_samples(phyloseq::sample_sums(ps_s) > 0, ps_s)
  if (phyloseq::nsamples(ps_s) < 1 || phyloseq::ntaxa(ps_s) < 2) {
    message("Skipping site ", s, " (not enough data)."); next
  }
  
  M <- as(phyloseq::otu_table(ps_s), "matrix")
  if (!phyloseq::taxa_are_rows(ps_s)) M <- t(M)
  
  keep <- order(rowSums(M, na.rm = TRUE), decreasing = TRUE)[seq_len(min(30, nrow(M)))]
  M <- M[keep, , drop = FALSE]
  rownames(M) <- label_for_taxa(ps_s, rownames(M))
  
  # Optional sample annotations if columns exist
  ann <- as(phyloseq::sample_data(ps_s), "data.frame")
  ann$SampleID <- rownames(ann)
  ann_df <- ann %>%
    dplyr::transmute(
      SampleID,
      Zone = if ("zone" %in% names(ann)) ann$zone else NA,
      Organic = if ("organic_pct" %in% names(ann)) suppressWarnings(as.numeric(ann$organic_pct)) else NA_real_,
      Survival = if ("survival_pct" %in% names(ann)) suppressWarnings(as.numeric(ann$survival_pct)) else NA_real_
    ) %>% tibble::column_to_rownames("SampleID")
  
  # align columns to annotation
  common_cols <- intersect(colnames(M), rownames(ann_df))
  M <- M[, common_cols, drop = FALSE]
  ann_df <- ann_df[common_cols, , drop = FALSE]
  
  out_png <- file.path("exports/heatmaps_by_site", paste0("heatmap_top30_", s, ".png"))
  pheatmap::pheatmap(
    M,
    scale = "row",
    clustering_distance_cols = "euclidean",
    clustering_method = "average",
    annotation_col = ann_df,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    main = paste0("Top 30 taxa — site ", s),
    filename = out_png,
    width = 11, height = 8
  )
}

# ---------- E) (Optional) per-site tree for those top 30 ----------
for (s in sites) {
  ps_s <- phyloseq::subset_samples(ps_tx_rel, site_abbrev == s)
  ps_s <- phyloseq::prune_samples(phyloseq::sample_sums(ps_s) > 0, ps_s)
  if (phyloseq::nsamples(ps_s) < 1 || phyloseq::ntaxa(ps_s) < 2) next
  
  M <- as(phyloseq::otu_table(ps_s), "matrix")
  if (!phyloseq::taxa_are_rows(ps_s)) M <- t(M)
  keep <- order(rowSums(M, na.rm = TRUE), decreasing = TRUE)[seq_len(min(30, nrow(M)))]
  ps_top <- phyloseq::prune_taxa(rownames(M)[keep], ps_s)
  
  out_tree <- file.path("exports/tree", paste0("tree_top30_", s, ".png"))
  png(out_tree, width = 1400, height = 1100, res = 200)
  plot(phyloseq::phy_tree(ps_top),
       main = paste0("Phylogenetic (proxy) tree — top 30 taxa at ", s))
  dev.off()
}

message("✅ Wrote: exports/tree/*.png and exports/heatmaps_by_site/heatmap_top30_*.png")

######
######
######
######

## ===================== PATCH: ORGANIC%, SURVIVAL%, ENV =====================

# Helper: robust numeric parser (handles numeric, factor, and "12.3%" strings)
to_num <- function(x){
  if (is.numeric(x)) return(as.numeric(x))
  if (is.factor(x))  x <- as.character(x)
  readr::parse_number(as.character(x))
}

# Helper: pick first existing column (define if not already present)
if (!exists("first_existing_col")) {
  first_existing_col <- function(df, cand){
    cc <- intersect(cand, names(df)); if (length(cc)) cc[1] else NULL
  }
}

# --- ORGANIC % from CHEM (e.g., "Percent Org Total") ---
# After clean_names(), likely becomes "percent_org_total"
ORGANIC_COL <- first_existing_col(
  chem,
  c("percent_org_total","percent_organic_total","percent_organic","percent_org","organic_percent")
)

# Create chem$organic_pct as numeric percent (3.21% -> 3.21)
chem$organic_pct <- if (!is.null(ORGANIC_COL)) to_num(chem[[ORGANIC_COL]]) else NA_real_

# Normalize key for the join
chem$site_abbrev <- toupper(trimws(as.character(chem$site_abbrev)))
mdf$site_abbrev  <- toupper(trimws(as.character(mdf$site_abbrev)))

# Merge organic % into MDF (by site_abbrev)
mdf <- dplyr::left_join(
  mdf,
  chem |>
    dplyr::select(site_abbrev, organic_pct,
                  cu = !!first_existing_col(chem, c("cu")), 
                  ni = !!first_existing_col(chem, c("ni")),
                  zn = !!first_existing_col(chem, c("zn")),
                  co = !!first_existing_col(chem, c("co"))),
  by = "site_abbrev"
)

# --- SURVIVAL % from MDF (e.g., "survival_100_percent_96hrs") ---
SURV_COL <- first_existing_col(
  mdf,
  c("survival_100_percent_96hrs","survival_percent","percent_survival","surv_pct","survival")
)

mdf$survival_pct <- if (!is.null(SURV_COL)) to_num(mdf[[SURV_COL]]) else NA_real_

# If survival is 0–1, convert to 0–100
if (is.finite(suppressWarnings(max(mdf$survival_pct, na.rm = TRUE))) &&
    suppressWarnings(max(mdf$survival_pct, na.rm = TRUE)) <= 1) {
  mdf$survival_pct <- mdf$survival_pct * 100
}

# --- Build env_keep safely (no rownames() inside dplyr) ---
mdf$SampleID <- if ("SampleID" %in% names(mdf) && !all(is.na(mdf$SampleID))) {
  as.character(mdf$SampleID)
} else {
  # fallback to rownames if needed
  if (!is.null(rownames(mdf))) rownames(mdf) else seq_len(nrow(mdf)) |> as.character()
}

env_keep <- mdf |>
  dplyr::transmute(
    SampleID,
    cu = suppressWarnings(as.numeric(cu)),
    ni = suppressWarnings(as.numeric(ni)),
    zn = suppressWarnings(as.numeric(zn)),
    co = suppressWarnings(as.numeric(co)),
    organic_pct = suppressWarnings(as.numeric(organic_pct)),
    survival_pct = suppressWarnings(as.numeric(survival_pct)),
    zone = droplevels(factor(zone)),
    site_abbrev = droplevels(factor(site_abbrev))
  )

# Diagnostics
message("✅ organic_pct non-NA: ", sum(is.finite(env_keep$organic_pct)))
message("✅ survival_pct non-NA: ", sum(is.finite(env_keep$survival_pct)))
readr::write_csv(env_keep, "exports/env_keep_for_CCA.csv")

# =================== END PATCH (continue with CCA/plots) ====================

# ===============================================================
# A) Detect and merge ORGANIC % and SURVIVAL % into metadata (mdf)
# ===============================================================

# --- Detect columns ---
ORGANIC_COL <- first_existing_col(chem, c("percent_org_total", "percent_organic", "percent_org", "percent_organics"))
SURV_COL    <- first_existing_col(mdf, c("survival_100_percent_96hrs", "survival_percent", "survival", "surv_96h"))

message("Detected organic column in chem: ", ORGANIC_COL %||% "None")
message("Detected survival column in sample data: ", SURV_COL %||% "None")

# --- Clean organic % from chem and merge into mdf ---
if(!is.null(ORGANIC_COL) && ORGANIC_COL %in% names(chem)){
  chem <- chem %>%
    mutate(
      !!ORGANIC_COL := readr::parse_number(!!sym(ORGANIC_COL))  # clean % signs to numeric
    )
  org_df <- chem %>% dplyr::select(site_abbrev, !!ORGANIC_COL)
  names(org_df)[2] <- "organic_pct"
  mdf <- mdf %>% dplyr::left_join(org_df, by = "site_abbrev")
} else {
  mdf$organic_pct <- NA_real_
}

# --- Clean survival % from sample data ---
mdf$survival_pct <- if(!is.null(SURV_COL) && SURV_COL %in% names(mdf)){
  readr::parse_number(mdf[[SURV_COL]])
} else NA_real_

# ===============================================================
# B) Build safe environmental matrix for CCA and ordinations
# ===============================================================

env_keep <- mdf %>%
  dplyr::transmute(
    SampleID = rownames(.),
    cu = as.numeric(cu),
    ni = as.numeric(ni),
    zn = as.numeric(zn),
    co = as.numeric(co),
    organic_pct = as.numeric(organic_pct),
    survival_pct = as.numeric(survival_pct),
    zone = droplevels(factor(zone)),
    site_abbrev = droplevels(factor(site_abbrev))
  )

# ===============================================================
# C) Confirm completeness
# ===============================================================
summary(env_keep)
readr::write_csv(env_keep, "exports/env_keep_for_CCA.csv")
message("✅ Environmental metadata ready for CCA and vector-based analyses.")

# ===============================================================
# D) Canonical Correspondence Analysis (CCA)
# ===============================================================

# Select taxa table (relative abundance)
ps_rel <- phyloseq::transform_sample_counts(ps, function(x) if(sum(x)==0) x else x/sum(x))
otu_rel <- as.data.frame(phyloseq::otu_table(ps_rel))
otu_rel <- otu_rel[, match(env_keep$SampleID, colnames(otu_rel)), drop = FALSE]

# Run CCA using contaminants + organics + survival
cca_model <- vegan::cca(t(otu_rel) ~ cu + ni + zn + co + organic_pct + survival_pct, data = env_keep)

# Plot CCA biplot
png("exports/cca_contaminants_organics_survival.png", width=1200, height=900, res=200)
plot(cca_model, display=c("sites","species","bp"), main="CCA — Contaminants, Organics & Survival")
dev.off()

# Export scores
scores_sites <- as.data.frame(scores(cca_model, display="sites"))
scores_species <- as.data.frame(scores(cca_model, display="species"))
readr::write_csv(scores_sites, "exports/cca_scores_sites.csv")
readr::write_csv(scores_species, "exports/cca_scores_species.csv")

# ===============================================================
# E) Phylogenetic Tree by Sample
# ===============================================================

tree <- phyloseq::phy_tree(ps)
if(!is.null(tree)){
  png("exports/phylogenetic_tree_by_sample.png", width=1400, height=1000, res=200)
  plot(tree, main="Phylogenetic Tree of Taxa per Sample")
  dev.off()
} else {
  message("⚠️ No phylogenetic tree in phyloseq object.")
}

# ===============================================================
# F) Heatmap: Survival × Contaminants × Dominant Taxa
# ===============================================================

# Identify most prevalent taxa
taxa_abund <- data.frame(taxa = taxa_names(ps), abundance = taxa_sums(ps))
top_taxa <- taxa_abund %>% arrange(desc(abundance)) %>% slice_head(n = 40) %>% pull(taxa)

# Create matrix of top taxa (relative abundance)
ps_top <- prune_taxa(top_taxa, ps_rel)
mat <- as.data.frame(otu_table(ps_top))
mat <- mat %>% mutate(Taxon = rownames(.)) %>% relocate(Taxon)

# Combine with survival and contaminants
meta_for_heat <- env_keep %>%
  dplyr::select(SampleID, cu, ni, zn, co, organic_pct, survival_pct)
rownames(meta_for_heat) <- meta_for_heat$SampleID

# Order samples by survival
mat_order <- mat[, c("Taxon", meta_for_heat$SampleID[order(meta_for_heat$survival_pct, decreasing = TRUE)])]

pheatmap::pheatmap(
  mat_order %>% column_to_rownames("Taxon"),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  annotation_col = meta_for_heat %>% select(survival_pct, organic_pct),
  main = "Heatmap: Top Taxa vs Contaminants & Survival"
)

# ===============================================================
# G) Alpha & Beta Diversity by Location
# ===============================================================

# Alpha diversity boxplots
alpha_df <- readr::read_csv("exports/alpha_metrics_merged.csv", show_col_types = FALSE)
p_alpha_loc <- ggplot(alpha_df, aes(x = site_abbrev, y = Shannon, fill = site_abbrev)) +
  geom_boxplot() + geom_jitter(width = 0.2) +
  labs(title = "Alpha Diversity (Shannon) by Location", x = "Location", y = "Shannon Index") +
  theme_bw(base_size = 12) + theme(legend.position="none")
ggsave("exports/alpha_diversity_by_location.png", p_alpha_loc, width=8, height=5, dpi=300)

# Beta diversity PCoA colored by SampleID
ord_pcoa <- ordinate(ps_rel, method="PCoA", distance="bray")
scores <- as.data.frame(ord_pcoa$vectors) %>% tibble::rownames_to_column("SampleID")
p_beta <- ggplot(scores, aes(Axis.1, Axis.2, color = SampleID)) +
  geom_point(size=2.2, alpha=0.8) +
  ggrepel::geom_text_repel(aes(label=SampleID), size=2.5, show.legend=FALSE) +
  theme_bw(base_size=12) + theme(legend.position="none") +
  labs(title="Beta Diversity (PCoA — Bray-Curtis)", x="PCoA1", y="PCoA2")
ggsave("exports/beta_diversity_by_sampleid.png", p_beta, width=8, height=6, dpi=300)

# ===============================================================
# ✅ END — Ready for integration into your main pipeline
# ===============================================================



######
######

# ---------- A) Identify organic % and survival columns ----------
# From chemistry: try to find an "organic %" proxy (e.g., LOI/TOC/etc.)
find_org_col <- function(df){
  cand <- names(df)
  # prioritize % or fractional organics names
  keys <- c("organic", "org", "loi", "toc", "t.o.c", "total_organic", "percent_organic", "%_organic")
  hits <- cand[Reduce(`|`, lapply(keys, function(k) grepl(k, cand, ignore.case = TRUE)))]
  # avoid false positives like "organization"
  hits <- hits[!grepl("organizat|organis", hits, ignore.case = TRUE)]
  hits[1]
}

# From metadata: try to find survival %
find_surv_col <- function(df){
  cand <- names(df)
  keys <- c("survival", "surv", "mortality", "percent_survival", "surv_pct", "%survival")
  hits <- cand[Reduce(`|`, lapply(keys, function(k) grepl(k, cand, ignore.case = TRUE)))]
  hits[1]
}

ORGANIC_COL <- find_org_col(chem)
SURV_COL    <- find_surv_col(mdf)

message("Detected organic column in chem: ", ORGANIC_COL %||% "None")
message("Detected survival column in sample data: ", SURV_COL %||% "None")

# Coerce numeric copies on mdf for modeling/annotations
mdf$organic_pct <- if (!is.null(ORGANIC_COL) && ORGANIC_COL %in% names(chem)) {
  # merge chem organic into mdf by site_abbrev
  org_df <- chem |> dplyr::select(site_abbrev, !!ORGANIC_COL)
  names(org_df)[2] <- "organic_raw"
  tmp <- mdf |> dplyr::left_join(org_df, by = "site_abbrev")
  suppressWarnings(as.numeric(tmp$organic_raw))
} else NA_real_

mdf$survival_pct <- if (!is.null(SURV_COL) && SURV_COL %in% names(mdf)) {
  suppressWarnings(as.numeric(mdf[[SURV_COL]]))
} else NA_real_

# A safe “env” frame we’ll reuse
env_keep <- mdf |>
  dplyr::transmute(
    SampleID = rownames(.),
    cu = as.numeric(cu), ni = as.numeric(ni), zn = as.numeric(zn), co = as.numeric(co),
    organic_pct = as.numeric(organic_pct),
    survival_pct = as.numeric(survival_pct),
    zone = droplevels(factor(zone)),
    site_abbrev = droplevels(factor(site_abbrev))
  )




# ====
# 

# ---------- B) CCA with arrows & labeled samples ----------
# Hellinger-transform the community (recommended for CCA)
ps_rel0 <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
comm    <- as(otu_table(ps_rel0), "matrix"); if(!taxa_are_rows(ps_rel0)) comm <- t(comm)
comm_hel <- vegan::decostand(comm, method = "hellinger")

# Align env to samples in comm
env2 <- env_keep |> dplyr::filter(SampleID %in% colnames(comm_hel))
comm_hel <- comm_hel[, env2$SampleID, drop = FALSE]

# Select numeric predictors that exist
vars_num <- c("cu","ni","zn","co","organic_pct","survival_pct")
X <- env2[, vars_num, drop = FALSE] |> dplyr::mutate(dplyr::across(everything(), as.numeric))

# Remove columns with all NA or zero variance
good_cols <- names(X)[sapply(X, function(v) sum(is.finite(v)) >= 3 && stats::var(v, na.rm=TRUE) > 0)]
X <- X[, good_cols, drop = FALSE]

stopifnot(nrow(env2) >= 3, ncol(comm_hel) >= 3, ncol(X) >= 1)

cca_fit <- vegan::cca(t(comm_hel) ~ ., data = X)  # sites = samples
# Save a quick summary
sink("exports/cca_summary.txt"); print(summary(cca_fit)); sink()

# Tidy site scores for plotting & labels
sc_sites <- as.data.frame(vegan::scores(cca_fit, display = "sites"))
sc_sites$SampleID <- rownames(sc_sites)
sc_sites <- sc_sites |>
  dplyr::left_join(mdf |> dplyr::select(SampleID = rownames(mdf), zone, site_abbrev), by = "SampleID")

# Tidy biplot arrows (env loadings)
sc_arrows <- as.data.frame(vegan::scores(cca_fit, display = "bp"))
sc_arrows$var <- rownames(sc_arrows)

library(ggrepel)
p_cca <- ggplot(sc_sites, aes(CCA1, CCA2, color = zone)) +
  geom_point(size = 2.4, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = 80, show.legend = FALSE) +
  geom_segment(data = sc_arrows,
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.6) +
  geom_text(data = sc_arrows, aes(CCA1, CCA2, label = var),
            vjust = -0.4, size = 3) +
  labs(title = "CCA — community vs contaminants, organic %, survival",
       color = "Zone") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave("exports/CCA_biplot_samples_labeled.png", p_cca, width = 9, height = 7, dpi = 300)

# 
# 
# 
## ---------- C) Tree of taxa with per-sample abundance heat (one panel per sample) ----------
# If a true phylogeny exists in `ps`, use it; else build a taxonomy-based proxy tree.
ps_tree <- ps
if (is.null(phyloseq::phy_tree(ps_tree, errorIfNULL = FALSE))) {
  message("No phylogenetic tree present; constructing a taxonomy-based proxy tree.")
  tt <- as.data.frame(phyloseq::tax_table(ps_tree))
  # Build a simple hierarchical string and cluster
  tax_path <- apply(tt, 1, function(x) paste(na.omit(as.character(x)), collapse = ";"))
  D <- stats::dist(stringdist::stringdistmatrix(tax_path, tax_path, method = "lv"))
  hc <- stats::hclust(D, method = "average")
  tree_proxy <- ape::as.phylo(hc)
  phyloseq::phy_tree(ps_tree) <- tree_proxy
}

# Make relative abundances and cap small values
ps_tree_rel <- phyloseq::transform_sample_counts(ps_tree, function(x) x/sum(x))
otu_rel <- as(otu_table(ps_tree_rel), "matrix"); if(!taxa_are_rows(ps_tree_rel)) otu_rel <- t(otu_rel)
otu_rel <- pmin(otu_rel, 0.2)

# Plot a multi-sample tree heat (one big heat), then a faceted PDF for per-sample
dir.create("exports/tree", showWarnings = FALSE)

# Overall tree with a compact heatmap strip per sample
pheatmap::pheatmap(
  otu_rel,
  show_rownames = FALSE, show_colnames = TRUE,
  filename = "exports/tree/taxa_tree_heat_overall.png", width = 14, height = 10
)

# Optional: one PNG per sample highlighting that sample’s abundance
samples_list <- colnames(otu_rel)
for (s in samples_list) {
  mat <- cbind(Abundance = otu_rel[, s, drop = FALSE])
  pheatmap::pheatmap(
    mat, show_rownames = FALSE, show_colnames = TRUE,
    filename = file.path("exports/tree", paste0("taxa_tree_heat_", s, ".png")),
    width = 5, height = 8
  )
}



###
###
###
# ---------- D) Heatmap of survival, contaminants, and top taxa ----------
# Get relative abundance by Genus (collapse safely)
rk <- intersect(c("Genus","genus","Family","family"), colnames(as.data.frame(tax_table(ps))))
use_rank <- if (length(rk)) rk[1] else NA_character_
ps_gen <- if (!is.na(use_rank)) suppressWarnings(tax_glom(ps, taxrank = use_rank, NArm = FALSE)) else ps
ps_gen_rel <- transform_sample_counts(ps_gen, function(x) x/sum(x))
otu_g <- as.data.frame(otu_table(ps_gen_rel)); if(!taxa_are_rows(ps_gen_rel)) otu_g <- t(otu_g)

# Pick top taxa by overall prevalence (sum of rel. abund)
topN <- 25
keep <- order(rowSums(otu_g, na.rm = TRUE), decreasing = TRUE)[seq_len(min(topN, nrow(otu_g)))]
mat_top <- as.matrix(otu_g[keep, , drop = FALSE])

# Taxon labels
tx <- as.data.frame(tax_table(ps_gen_rel))
lab <- if (!is.na(use_rank) && use_rank %in% names(tx)) tx[rownames(mat_top), use_rank] else rownames(mat_top)
rownames(mat_top) <- make.names(lab, unique = TRUE)

# Annotations for samples
ann <- mdf |>
  dplyr::transmute(
    SampleID = rownames(.),
    Survival = survival_pct,
    Cu = as.numeric(cu), Ni = as.numeric(ni), Zn = as.numeric(zn), Co = as.numeric(co),
    Organic = as.numeric(organic_pct),
    Zone = zone, Location = site_abbrev
  ) |>
  tibble::column_to_rownames("SampleID")

mat_top <- mat_top[, rownames(ann), drop = FALSE]  # ensure same order

pheatmap::pheatmap(
  mat_top,
  annotation_col = ann[, c("Survival","Cu","Ni","Zn","Co","Organic","Zone","Location")],
  scale = "row",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  show_rownames = TRUE,
  filename = "exports/heatmap_top_taxa_survival_contaminants.png",
  width = 14, height = 9
)


# 

# ---------- E) Alpha by location (site_abbrev) ----------
alpha_df <- readr::read_csv("exports/alpha_metrics_merged.csv", show_col_types = FALSE)
alpha_df$site_abbrev <- factor(alpha_df$site_abbrev)

p_alpha_loc <- ggplot(alpha_df, aes(site_abbrev, Shannon, fill = site_abbrev)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  labs(x = "Location (site_abbrev)", y = "Shannon",
       title = "Alpha diversity (Shannon) by location") +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  coord_flip()
ggsave("exports/alpha_Shannon_by_location.png", p_alpha_loc, width = 7, height = 6, dpi = 300)

# 
# ---------- F) PCoA with sample labels ----------
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))
dist_bray <- phyloseq::distance(ps_rel, method = "bray")
ord      <- ordinate(ps_rel, method = "PCoA", distance = dist_bray)

scores <- as.data.frame(ord$vectors)
scores$SampleID <- rownames(scores)
scores <- scores |>
  dplyr::left_join(mdf |> dplyr::select(SampleID = rownames(mdf), zone, site_abbrev), by = "SampleID")

pc1 <- round(ord$values$Relative_eig[1]*100,1)
pc2 <- round(ord$values$Relative_eig[2]*100,1)

p_lab <- ggplot(scores, aes(Axis.1, Axis.2, color = zone, label = SampleID)) +
  geom_point(size = 2.4, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 80, show.legend = FALSE) +
  labs(x = glue::glue("PCoA1 ({pc1}%)"), y = glue::glue("PCoA2 ({pc2}%)"),
       title = "PCoA (Bray-Curtis) — labeled by sample ID") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave("exports/beta_pcoa_labeled_by_sample.png", p_lab, width = 9, height = 7, dpi = 300)



###
###
###
# ---------- G) PCA on top taxa with envfit arrows ----------
# Use the same top taxa matrix from D (mat_top), apply Hellinger transform again just in case
mat_pca <- vegan::decostand(t(mat_top), method = "hellinger")  # samples x taxa
pca_fit <- stats::prcomp(mat_pca, center = TRUE, scale. = FALSE)

# Site scores
pca_sites <- as.data.frame(pca_fit$x) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::left_join(env_keep, by = "SampleID")

# Variance %
varp <- round(100 * summary(pca_fit)$importance["Proportion of Variance", 1:2], 1)

# Fit env arrows (contaminants + organic)
env_for_fit <- pca_sites |> dplyr::select(cu,ni,zn,co,organic_pct)
ef <- vegan::envfit(as.matrix(pca_fit$x[,1:2]), env_for_fit, permutations = 9999)

ef_ar <- as.data.frame(ef$vectors$arrows * rep(ef$vectors$r, each = 2))
names(ef_ar) <- c("PC1","PC2")
ef_ar$var <- rownames(ef$vectors$arrows)

p_taxa_pca <- ggplot(pca_sites, aes(PC1, PC2, color = zone)) +
  geom_point(size = 2.4, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = 80, show.legend = FALSE) +
  geom_segment(data = ef_ar,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.6) +
  geom_text(data = ef_ar, aes(PC1, PC2, label = var), vjust = -0.4, size = 3) +
  labs(title = "PCA of most prevalent taxa (+ envfit: contaminants & organic)",
       x = glue::glue("PC1 ({varp[1]}%)"), y = glue::glue("PC2 ({varp[2]}%)"),
       color = "Zone") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave("exports/pca_top_taxa_envfit.png", p_taxa_pca, width = 9, height = 7, dpi = 300)

# ===========================================================


nrow(readr::read_csv("seqtab_ci.csv", show_col_types = FALSE))
nrow(readxl::read_excel(meta_path))


