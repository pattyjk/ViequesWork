############################################################
# ITS: Contaminants + Taxa Relationship Plots (ONE SCRIPT)
# Labeled by Site (Site abbreviation), Colored by Zone
# Exports figures + tables to exports_contam_taxa/
############################################################

# ---- Packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggrepel)
  library(indicspecies)
  library(permute)
})

# Optional (nice heatmaps). If not installed, script will still run other plots.
has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)

# -----------------------------
# 0) Paths
# -----------------------------
wd <- "C:/Users/valer/OneDrive/Desktop/Amplicon_P.R_Sediments/ITS_R/ITS_ampl_sed"
setwd(wd)

out_dir <- file.path(getwd(), "exports_contam_taxa")
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE)

# -----------------------------
# 1) Load data
# -----------------------------
meta <- read.csv("cleaned_metadata.csv", check.names = FALSE)

stopifnot(file.exists("asv.table.txt"))
asv_raw <- read.delim("asv.table.txt", check.names = FALSE)

# Community matrix: rows = samples, cols = ASVs (your file is 29 x 2380 already)
asv <- as.matrix(asv_raw)
mode(asv) <- "numeric"

# -----------------------------
# 2) Helpers
# -----------------------------
clean_id <- function(x){
  x %>%
    as.character() %>%
    stringr::str_replace_all("\u00A0"," ") %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("\\s+","") %>%
    toupper()
}

num_clean <- function(x){
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "%", "")
  x <- stringr::str_replace_all(x, ",", "")
  x <- stringr::str_replace_all(x, "<", "")
  x <- stringr::str_replace_all(x, "ND|N/D|n/d", NA_character_)
  x <- stringr::str_replace_all(x, "—|–|NA|NaN|Inf|-Inf", NA_character_)
  x[x == ""] <- NA_character_
  suppressWarnings(as.numeric(x))
}

# ---- Align sample IDs between asv and metadata
# asv rownames come from file rownames; ensure they match meta SampleID
rownames(asv) <- clean_id(rownames(asv))
meta <- meta %>% mutate(SampleID_clean = clean_id(SampleID))

common <- intersect(rownames(asv), meta$SampleID_clean)
if(length(common) < 5) {
  cat("ASV rownames head:", paste(head(rownames(asv), 10), collapse=", "), "\n")
  cat("Meta SampleID head:", paste(head(meta$SampleID_clean, 10), collapse=", "), "\n")
  stop("Too few overlapping SampleIDs between asv.table.txt and metadata.")
}

asv  <- asv[common, , drop=FALSE]
meta2 <- meta %>%
  filter(SampleID_clean %in% common) %>%
  arrange(match(SampleID_clean, rownames(asv)))

# ---- Variables (your column names)
meta2 <- meta2 %>%
  mutate(
    Zone = as.factor(Zone),
    Location = as.factor(Location),
    SiteAbbrev = as.factor(`Site  abbreviation`)
  )

# -----------------------------
# 3) Zone color template (EDIT HERE if you have specific hex codes)
# -----------------------------
# If you have 4 zones, set exact hex codes here in the order of levels(meta2$Zone)
# Example placeholders:
# Ensure Zone is a factor
meta2$Zone <- as.factor(meta2$Zone)
zone_colors <- zone_colors[levels(meta2$Zone)]
# Your custom palette
zone_colors <- c(
  "1" = "#C6A500",  # mustard
  "2" = "#2E8B57",  # green
  "3" = "#E67E22",  # orange
  "4" = "#C0392B"   # red
)

# Check alignment
print(levels(meta2$Zone))
print(zone_colors)
# -----------------------------
# 4) Environmental variables (two models)
#    Model A: n=29 exclude Ni
#    Model B: n=26 include Ni (drops LAN1-3 if Ni missing)
# -----------------------------
env <- meta2 %>%
  select(`Percent Org Total`, `96hrs Survial %`, `Cu Average`, `Ni Average`, `Zn Average`, `Co Average`) %>%
  mutate(across(everything(), num_clean))
rownames(env) <- rownames(asv)

bad_env <- rownames(env)[!complete.cases(env)]
write.csv(env[bad_env, , drop=FALSE], file.path(out_dir,"tables","env_missing_rows.csv"))

env_noNi <- env %>% select(-`Ni Average`)
env_noNi_scaled <- as.data.frame(scale(env_noNi))

keep_withNi <- complete.cases(env)
env_withNi <- env[keep_withNi, , drop=FALSE]
env_withNi_scaled <- as.data.frame(scale(env_withNi))

# -----------------------------
# 5) Community transforms + distance (use filtering to avoid “stuck” ordinations)
# -----------------------------
# Filter rare ASVs to reduce sparsity for ordination/IndVal
# (present in >= 3 samples; adjust to 2 if you want less strict)
prev <- colSums(asv > 0)
asv_f <- asv[, prev >= 3, drop=FALSE]
cat("ASVs kept after prevalence filter:", ncol(asv_f), "of", ncol(asv), "\n")
write.csv(
  data.frame(ASV_Sequence = colnames(asv), Prevalence = prev),
  file.path(out_dir,"tables","ASV_prevalence.csv"),
  row.names = FALSE
)

# Relative abundance + sqrt transform stabilizes Bray for sparse ITS
rel_f  <- vegan::decostand(asv_f, "total")
rel_sqrt <- sqrt(rel_f)

# For RDA: Hellinger on filtered matrix
hell_f <- vegan::decostand(asv_f, "hellinger")

# -----------------------------
# 6) NMDS plots (labeled by SiteAbbrev, colored by Zone)
# -----------------------------
make_nmds_plot <- function(rel_mat, meta_df, env_scaled, label, suffix){
  bray <- vegan::vegdist(rel_mat, "bray")
  set.seed(123)
  nmds <- vegan::metaMDS(bray, k=2, trymax=200, autotransform=FALSE, trace=FALSE)
  
  scr <- as.data.frame(scores(nmds))
  scr$SampleID  <- rownames(scr)
  scr$Zone      <- meta_df$Zone[match(scr$SampleID, rownames(rel_mat))]
  scr$SiteAbbrev<- meta_df$SiteAbbrev[match(scr$SampleID, rownames(rel_mat))]
  
  p <- ggplot(scr, aes(NMDS1, NMDS2, color=Zone)) +
    geom_point(size=3) +
    ggrepel::geom_text_repel(aes(label=SiteAbbrev), size=3, max.overlaps=200) +
    scale_color_manual(values=zone_colors) +
    theme_minimal(base_size=14) +
    labs(title=paste0("ITS NMDS (Bray–Curtis) – ", label),
         subtitle="Labels = Site abbreviation; Colors = Zone")
  
  ggsave(file.path(out_dir,"figures", paste0("NMDS_", suffix, "_labelSiteAbbrev_colorZone.png")),
         p, width=9, height=7, dpi=300)
  
  # envfit vectors
  fit <- vegan::envfit(nmds, env_scaled, permutations=999)
  vec <- data.frame(
    Variable = rownames(fit$vectors$arrows),
    NMDS1 = fit$vectors$arrows[,1],
    NMDS2 = fit$vectors$arrows[,2],
    r2 = fit$vectors$r,
    p  = fit$vectors$pvals
  ) %>% arrange(p)
  
  write.csv(vec, file.path(out_dir,"tables", paste0("envfit_vectors_", suffix, ".csv")), row.names=FALSE)
  
  # Add vectors (only p <= 0.05) for clarity
  vec_sig <- vec %>% filter(p <= 0.05)
  if(nrow(vec_sig) > 0) {
    p_vec <- p +
      geom_segment(data=vec_sig, aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
                   inherit.aes=FALSE, arrow=arrow(length=unit(0.02,"npc")), linewidth=0.6) +
      geom_text(data=vec_sig, aes(x=NMDS1, y=NMDS2, label=Variable),
                inherit.aes=FALSE, size=4)
    ggsave(file.path(out_dir,"figures", paste0("NMDS_envfitVectors_", suffix, ".png")),
           p_vec, width=9, height=7, dpi=300)
  }
  
  return(list(nmds=nmds, vec=vec))
}

# Model A (n=29, no Ni): use all rows of rel_sqrt
res_noNi <- make_nmds_plot(rel_sqrt, meta2, env_noNi_scaled,
                           label="Model A (n=29; Ni excluded)", suffix="noNi_n29")

# Model B (n=26, with Ni): subset rel_sqrt to complete-case rows
rel_sqrt_withNi <- rel_sqrt[rownames(env_withNi_scaled), , drop=FALSE]
meta_withNi <- meta2 %>% filter(SampleID_clean %in% rownames(rel_sqrt_withNi)) %>%
  arrange(match(SampleID_clean, rownames(rel_sqrt_withNi)))

res_withNi <- make_nmds_plot(rel_sqrt_withNi, meta_withNi, env_withNi_scaled,
                             label="Model B (n=26; Ni included)", suffix="withNi_n26")

# -----------------------------
# 7) RDA plots (filtered ASVs) – save base R PNG
# -----------------------------
# RDA Model A (no Ni) uses n=29
rda_noNi <- vegan::rda(hell_f ~ ., data=env_noNi_scaled)
write.csv(as.data.frame(vegan::anova.cca(rda_noNi, permutations=999)),
          file.path(out_dir,"tables","RDA_noNi_global_n29.csv"))
write.csv(as.data.frame(vegan::anova.cca(rda_noNi, by="term", permutations=999)),
          file.path(out_dir,"tables","RDA_noNi_byTerm_n29.csv"))

png(file.path(out_dir,"figures","RDA_noNi_n29.png"), width=1200, height=900, res=150)
plot(rda_noNi, scaling=2, main="ITS RDA (Hellinger, filtered ASVs) – n=29 (Ni excluded)")
dev.off()

# RDA Model B (with Ni) uses n=26 subset
hell_withNi <- hell_f[rownames(env_withNi_scaled), , drop=FALSE]
rda_withNi <- vegan::rda(hell_withNi ~ ., data=env_withNi_scaled)
write.csv(as.data.frame(vegan::anova.cca(rda_withNi, permutations=999)),
          file.path(out_dir,"tables","RDA_withNi_global_n26.csv"))
write.csv(as.data.frame(vegan::anova.cca(rda_withNi, by="term", permutations=999)),
          file.path(out_dir,"tables","RDA_withNi_byTerm_n26.csv"))

png(file.path(out_dir,"figures","RDA_withNi_n26.png"), width=1200, height=900, res=150)
plot(rda_withNi, scaling=2, main="ITS RDA (Hellinger, filtered ASVs) – n=26 (Ni included)")
dev.off()

# -----------------------------
# 8) “Normal” contaminant plots (scatter + bubble), labeled by SiteAbbrev, colored by Zone
# -----------------------------
meta_plot <- meta2 %>%
  mutate(across(c(`Percent Org Total`,`96hrs Survial %`,`Cu Average`,`Ni Average`,`Zn Average`,`Co Average`), num_clean))

# Survival vs each metal
metal_cols <- c("Cu Average","Ni Average","Zn Average","Co Average")
for(m in metal_cols){
  p <- ggplot(meta_plot, aes(x=.data[[m]], y=`96hrs Survial %`, color=Zone, label=SiteAbbrev)) +
    geom_point(size=3) +
    ggrepel::geom_text_repel(size=3, max.overlaps=200) +
    scale_color_manual(values=zone_colors) +
    theme_minimal(base_size=14) +
    labs(title=paste0("Survival vs ", m),
         subtitle="Labels=Site abbreviation; Colors=Zone",
         x=m, y="96hr Survival (%)")
  ggsave(file.path(out_dir,"figures", paste0("Survival_vs_", gsub(" ","_",m), ".png")),
         p, width=9, height=7, dpi=300)
}

# Bubble plot: Cu vs Zn, size = Survival
p_bub <- ggplot(meta_plot, aes(`Cu Average`, `Zn Average`, size=`96hrs Survial %`, color=Zone, label=SiteAbbrev)) +
  geom_point(alpha=0.85) +
  ggrepel::geom_text_repel(size=3, max.overlaps=200) +
  scale_color_manual(values=zone_colors) +
  theme_minimal(base_size=14) +
  labs(title="Metal landscape: Cu vs Zn (bubble size = 96hr Survival)",
       subtitle="Labels=Site abbreviation; Colors=Zone",
       x="Cu Average", y="Zn Average", size="Survival (%)")
ggsave(file.path(out_dir,"figures","Bubble_Cu_vs_Zn_sizeSurvival.png"),
       p_bub, width=10, height=7, dpi=300)


#rda labeled 

library(ggrepel)

make_rda_plot <- function(rda_model, meta_df, env_scaled, label, suffix){
  
  # Extract site scores
  sites <- as.data.frame(scores(rda_model, display="sites", scaling=2))
  sites$SampleID <- rownames(sites)
  sites$Zone <- meta_df$Zone[match(sites$SampleID, meta_df$SampleID_clean)]
  sites$SiteAbbrev <- meta_df$SiteAbbrev[match(sites$SampleID, meta_df$SampleID_clean)]
  
  # Extract environmental vectors
  env_vec <- as.data.frame(scores(rda_model, display="bp", scaling=2))
  env_vec$Variable <- rownames(env_vec)
  
  # Variance explained
  var_expl <- summary(rda_model)$cont$importance[2,1:2] * 100
  xlab <- paste0("RDA1 (", round(var_expl[1],1), "%)")
  ylab <- paste0("RDA2 (", round(var_expl[2],1), "%)")
  
  p <- ggplot() +
    geom_point(data=sites,
               aes(RDA1, RDA2, color=Zone),
               size=3) +
    ggrepel::geom_text_repel(data=sites,
                             aes(RDA1, RDA2, label=SiteAbbrev),
                             size=3,
                             max.overlaps=200) +
    geom_segment(data=env_vec,
                 aes(x=0, y=0, xend=RDA1, yend=RDA2),
                 arrow=arrow(length=unit(0.02,"npc")),
                 linewidth=0.8) +
    geom_text(data=env_vec,
              aes(RDA1, RDA2, label=Variable),
              size=4,
              vjust=-0.5) +
    scale_color_manual(values=zone_colors) +
    theme_minimal(base_size=14) +
    labs(title=paste0("ITS RDA – ", label),
         subtitle="Points = Sites (Site Abbrev); Arrows = Environmental Variables",
         x=xlab,
         y=ylab)
  
  ggsave(file.path(out_dir,"figures",
                   paste0("RDA_", suffix, "_labelSiteAbbrev_colorZone.png")),
         p, width=9, height=7, dpi=300)
  
  return(p)
}

make_rda_plot(
  rda_withNi,
  meta_withNi,
  env_withNi_scaled,
  label = "Model B (n=26; Ni included)",
  suffix = "withNi_n26"
)
make_rda_plot(
  rda_noNi,
  meta2,
  env_noNi_scaled,
  label = "Model A (n=29; Ni excluded)",
  suffix = "noNi_n29"
)
#
#

make_rda_plot2 <- function(rda_model, meta_df, label, suffix,
                           arrow_mult = 0.6,   # shrink arrows (0.4–0.8 usually best)
                           point_size = 3,
                           site_label_size = 3,
                           env_label_size = 4){
  
  # Site scores
  sites <- as.data.frame(vegan::scores(rda_model, display="sites", scaling=2))
  sites$SampleID <- rownames(sites)
  
  # Match metadata (SampleID_clean must exist)
  sites$Zone <- meta_df$Zone[match(sites$SampleID, meta_df$SampleID_clean)]
  sites$SiteAbbrev <- meta_df$SiteAbbrev[match(sites$SampleID, meta_df$SampleID_clean)]
  
  # Env vectors (biplot scores)
  env_vec <- as.data.frame(vegan::scores(rda_model, display="bp", scaling=2))
  env_vec$Variable <- rownames(env_vec)
  
  # Rescale arrows so they fit inside the cloud
  env_vec <- env_vec %>%
    mutate(RDA1 = RDA1 * arrow_mult,
           RDA2 = RDA2 * arrow_mult)
  
  # % variance explained
  var_expl <- summary(rda_model)$cont$importance[2,1:2] * 100
  xlab <- paste0("RDA1 (", round(var_expl[1],1), "%)")
  ylab <- paste0("RDA2 (", round(var_expl[2],1), "%)")
  
  p <- ggplot() +
    geom_point(data=sites, aes(RDA1, RDA2, color=Zone), size=point_size) +
    
    # Site labels (repelled)
    ggrepel::geom_text_repel(
      data=sites,
      aes(RDA1, RDA2, label=SiteAbbrev, color=Zone),
      size=site_label_size,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0
    ) +
    
    # arrows
    geom_segment(
      data=env_vec,
      aes(x=0, y=0, xend=RDA1, yend=RDA2),
      arrow=arrow(length=unit(0.02,"npc")),
      linewidth=0.8,
      color="black"
    ) +
    
    # Env labels (repelled)
    ggrepel::geom_text_repel(
      data=env_vec,
      aes(RDA1, RDA2, label=Variable),
      size=env_label_size,
      color="black",
      box.padding = 0.4,
      point.padding = 0.3,
      min.segment.length = 0
    ) +
    
    scale_color_manual(values=zone_colors) +
    theme_minimal(base_size=14) +
    theme(
      legend.position="right",
      plot.title = element_text(face="bold"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title=paste0("ITS RDA – ", label),
      subtitle="Points = Sites (Site Abbrev); Arrows = Environmental Variables",
      x=xlab, y=ylab, color="Zone"
    )
  
  # Export PNG (taller = less crowding)
  ggsave(
    filename = file.path(out_dir, "figures", paste0("RDA_", suffix, "_labelSiteAbbrev_colorZone.png")),
    plot = p,
    width = 12, height = 6, dpi = 300
  )
  
  return(p)
}

p_rda_A <- make_rda_plot2(
  rda_noNi,
  meta2,
  label = "Model A (n=29; Ni excluded)",
  suffix = "noNi_n29",
  arrow_mult = 0.55
)
p_rda_B <- make_rda_plot2(
  rda_withNi,
  meta_withNi,
  label = "Model B (n=26; Ni included)",
  suffix = "withNi_n26",
  arrow_mult = 0.55
)


# -----------------------------
#
# ============================================================
# 9) INDICATORS (GENUS-LEVEL) — outputs are Genus names
# ============================================================

library(indicspecies)
library(permute)
library(dplyr)
library(tidyr)

# --- Make sure tax_mat exists
# tax_mat <- readRDS("ITS_tax_mat.rds")

# --- Build GENUS relative abundance matrix from rel_f (samples x ASVs)
# rel_f must already exist and be samples x ASVs (ASV IDs are sequences)
asv_ids <- colnames(rel_f)
tax_ids <- rownames(tax_mat)

keep_asvs <- intersect(asv_ids, tax_ids)
cat("ASVs in rel_f:", length(asv_ids), "\n")
cat("ASVs with taxonomy:", length(keep_asvs), "\n")

if(length(keep_asvs) < 10) stop("Too few ASVs matched to taxonomy. Check ASV IDs vs tax_mat rownames.")

rel_tax <- rel_f[, keep_asvs, drop=FALSE]
tax_sub <- tax_mat[keep_asvs, , drop=FALSE]

genus_vec <- as.character(tax_sub[,"Genus"])
genus_vec[is.na(genus_vec) | genus_vec=="" | genus_vec=="NA"] <- "Unassigned"

# Collapse ASVs -> Genus (sum rel. abundances within genus)
genus_rel <- rowsum(t(rel_tax), group = genus_vec)
genus_rel <- t(genus_rel)  # samples x genera

# Optional genus filter (recommended: keep genera present in >=2 samples)
gen_prev <- colSums(genus_rel > 0)
genus_rel <- genus_rel[, gen_prev >= 1 drop=FALSE]
cat("Genera kept (prev>=2):", ncol(genus_rel), "\n")

# Save prevalence table for reporting
gen_prev_tbl <- data.frame(
  Genus = names(gen_prev),
  Prevalence = as.integer(gen_prev),
  MeanRelAbund = colMeans(genus_rel, na.rm=TRUE)
) %>% arrange(desc(Prevalence), desc(MeanRelAbund))

write.csv(gen_prev_tbl, file.path(out_dir,"tables","Genus_prevalence_meanRelAbund.csv"), row.names=FALSE)

# --- Ensure sample alignment to metadata
stopifnot(all(meta2$SampleID_clean %in% rownames(genus_rel)))
genus_rel <- genus_rel[meta2$SampleID_clean, , drop=FALSE]

# --- Genus-level IndVal function
run_indval_genus <- function(comm_rel, group_vec, label){
  g <- as.factor(group_vec)
  tab <- table(g)
  
  if(length(tab) < 2) return(NULL)
  if(all(tab == 1)) {
    message("IndVal skipped for ", label, " (all singleton groups).")
    return(NULL)
  }
  
  set.seed(123)
  ind <- indicspecies::multipatt(
    comm_rel, g,
    func = "IndVal.g",
    control = permute::how(nperm = 999)
  )
  
  all_df <- as.data.frame(ind$sign)
  all_df$Genus <- rownames(all_df)
  
  # Ensure p.value exists
  if(!("p.value" %in% names(all_df))) {
    pc <- intersect(names(all_df), c("p.value","p","pvalue"))
    if(length(pc)>0) names(all_df)[names(all_df)==pc[1]] <- "p.value"
  }
  
  sig_df <- all_df %>%
    filter(!is.na(p.value), p.value < 0.1) %>%
    arrange(p.value) %>%
    mutate(Grouping = label)
  
  list(all = all_df, sig = sig_df)
}

# --- Run genus-level indicators
indG_zone <- run_indval_genus(genus_rel, meta2$Zone, "Zone")
indG_loc  <- run_indval_genus(genus_rel, meta2$Location, "Location")
indG_site <- run_indval_genus(genus_rel, meta2$SiteAbbrev, "SiteAbbrev")

# --- Export (GENUS names)
if(!is.null(indG_zone)){
  write.csv(indG_zone$all, file.path(out_dir,"tables","IndVal_GENUS_Zone_ALL.csv"), row.names=FALSE)
  write.csv(indG_zone$sig, file.path(out_dir,"tables","IndVal_GENUS_Zone_SIG.csv"), row.names=FALSE)
}
if(!is.null(indG_loc)){
  write.csv(indG_loc$all,  file.path(out_dir,"tables","IndVal_GENUS_Location_ALL.csv"), row.names=FALSE)
  write.csv(indG_loc$sig,  file.path(out_dir,"tables","IndVal_GENUS_Location_SIG.csv"), row.names=FALSE)
}
if(!is.null(indG_site)){
  write.csv(indG_site$all, file.path(out_dir,"tables","IndVal_GENUS_SiteAbbrev_ALL.csv"), row.names=FALSE)
  write.csv(indG_site$sig, file.path(out_dir,"tables","IndVal_GENUS_SiteAbbrev_SIG.csv"), row.names=FALSE)
}

# --- Plot: Top 20 Genus indicators by Zone
if(!is.null(indG_zone) && nrow(indG_zone$sig) > 0){
  topG <- indG_zone$sig %>% slice_min(p.value, n=20)
  
  p_indG <- ggplot(topG, aes(x=reorder(Genus, stat), y=stat)) +
    geom_col() +
    coord_flip() +
    theme_minimal(base_size=13) +
    labs(title="Top 20 Fungal Genus Indicators (IndVal.g) – Zone",
         x="Genus", y="IndVal.g statistic")
  
  ggsave(file.path(out_dir,"figures","IndVal_GENUS_Top20_Zone.png"),
         p_indG, width=10, height=7, dpi=300)
}

out_dir
list.files(file.path(out_dir, "figures"))
p_indG_
nrow(indG_zone$sig)

summary(indG_zone$all$p.value)
nrow(indG_zone$all)
nrow(indG_zone$sig)
# Ensure alignment
genus_rel <- genus_rel[meta2$SampleID_clean, , drop=FALSE]

indG_loc <- run_indval_genus(genus_rel, meta2$Location, "Location")

# Check results
summary(indG_loc$all$p.value)
nrow(indG_loc$sig)



#
# ============================================================
# A) Top 20 Genus Indicators (Location) — ALWAYS makes a plot
# ============================================================
# ============================
# Top 20 Genus Indicators – Location (fixed)
# ============================

# Export indval tables (Location)
write.csv(indG_loc$all,
          file.path(out_dir,"tables","IndVal_GENUS_Location_ALL.csv"),
          row.names=FALSE)

write.csv(indG_loc$sig,
          file.path(out_dir,"tables","IndVal_GENUS_Location_SIG.csv"),
          row.names=FALSE)

# Top 20 (or fewer if <20 rows)
top_n <- min(20, nrow(indG_loc$all))

top20_loc <- indG_loc$all %>%
  filter(!is.na(p.value)) %>%
  arrange(p.value) %>%
  slice_head(n = top_n) %>%
  mutate(Significant = ifelse(p.value < 0.05, "p < 0.05", "n.s."))

write.csv(top20_loc,
          file.path(out_dir,"tables","IndVal_GENUS_Location_Top20.csv"),
          row.names=FALSE)

p_top20_loc <- ggplot(top20_loc,
                      aes(x=reorder(Genus, stat),
                          y=stat,
                          fill=Significant)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c("p < 0.05"="#C0392B",
                             "n.s."="#95A5A6")) +
  theme_minimal(base_size=14) +
  labs(title="Top Fungal Genus Indicators (IndVal.g) – Location",
       subtitle="Red = significant at p < 0.05 (top 20 by smallest p)",
       x="Genus",
       y="IndVal.g statistic")

ggsave(file.path(out_dir,"figures","IndVal_GENUS_Top20_Location.png"),
       p_top20_loc, width=10, height=7, dpi=300)

p_top20_loc
# 
# # 
# ============================================================
# Genus vs Contaminants Correlation Heatmap
# ============================================================

envG <- env_noNi
rownames(envG) <- meta2$SampleID_clean

# Align
genus_corr <- genus_rel[rownames(envG), , drop=FALSE]

# Keep top 30 genera by mean abundance
topG <- 30
top_genera <- names(sort(colMeans(genus_corr),
                         decreasing=TRUE))[1:min(topG, ncol(genus_corr))]

genus_top <- genus_corr[, top_genera, drop=FALSE]

# Spearman correlation
cor_mat_genus <- sapply(colnames(envG), function(v){
  apply(genus_top, 2, function(x){
    suppressWarnings(cor(x, envG[[v]],
                         method="spearman",
                         use="pairwise.complete.obs"))
  })
})
cor_mat_genus <- t(cor_mat_genus)

write.csv(cor_mat_genus,
          file.path(out_dir,"tables","SpearmanCorr_env_vs_topGenera.csv"))

if(requireNamespace("pheatmap", quietly = TRUE)){
  pheatmap::pheatmap(
    cor_mat_genus,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#2C3E50","white","#C0392B"))(50),
    main = "Spearman Correlation: Contaminants vs Top Fungal Genera",
    filename = file.path(out_dir,"figures","Heatmap_Spearman_env_vs_topGenera.png"),
    width = 10,
    height = 6
  )
}
# 
# 
# ============================================================
# B) ASV vs contaminants correlation heatmap (Spearman)
# ============================================================

# Choose env set:
# Use env_noNi for n=29 (keeps all samples), or env_withNi for n=26
envA <- env_noNi   # <- recommended for full n=29
rownames(envA) <- rownames(asv)  # match sample IDs if needed

# Use filtered rel abund (rel_f) for ASV correlations
rel_corr <- rel_f[rownames(envA), , drop=FALSE]

# Keep top 100 abundant ASVs (adjust to 50/200 if you want)
topN <- 100
top_asvs <- names(sort(colSums(rel_corr), decreasing=TRUE))[1:min(topN, ncol(rel_corr))]
rel_top <- rel_corr[, top_asvs, drop=FALSE]

# Spearman correlations: env vars (rows) x ASVs (cols)
cor_mat_asv <- sapply(colnames(envA), function(v){
  apply(rel_top, 2, function(x){
    suppressWarnings(cor(x, envA[[v]], method="spearman", use="pairwise.complete.obs"))
  })
})
cor_mat_asv <- t(cor_mat_asv)  # env vars as rows

write.csv(cor_mat_asv, file.path(out_dir,"tables","SpearmanCorr_env_vs_topASVs.csv"))

# Heatmap (requires pheatmap). If not installed, install.packages("pheatmap")
if(requireNamespace("pheatmap", quietly = TRUE)){
  pheatmap::pheatmap(
    cor_mat_asv,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = paste0("Spearman: Env vars vs Top ", topN, " ASVs (ITS)"),
    filename = file.path(out_dir,"figures","Heatmap_Spearman_env_vs_topASVs.png"),
    width = 12, height = 6
  )
} else {
  message("pheatmap not installed. CSV saved; install.packages('pheatmap') for PNG heatmap.")
}
#
# ============================================================
# C) Genus vs contaminants correlation heatmap (Spearman)
# ============================================================

# Use the same env set as above (n=29 recommended)
envG <- env_noNi
rownames(envG) <- meta2$SampleID_clean

# Align rows
genus_corr <- genus_rel[rownames(envG), , drop=FALSE]

# Keep top 30 genera by mean relative abundance (change to 50 if you want)
topG <- 30
top_genera <- names(sort(colMeans(genus_corr), decreasing=TRUE))[1:min(topG, ncol(genus_corr))]
genus_top <- genus_corr[, top_genera, drop=FALSE]

# Correlations: env vars (rows) x genera (cols)
cor_mat_genus <- sapply(colnames(envG), function(v){
  apply(genus_top, 2, function(x){
    suppressWarnings(cor(x, envG[[v]], method="spearman", use="pairwise.complete.obs"))
  })
})
cor_mat_genus <- t(cor_mat_genus)

write.csv(cor_mat_genus, file.path(out_dir,"tables","SpearmanCorr_env_vs_topGenera.csv"))

if(requireNamespace("pheatmap", quietly = TRUE)){
  pheatmap::pheatmap(
    cor_mat_genus,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = paste0("Spearman: Env vars vs Top ", topG, " Genera (ITS)"),
    filename = file.path(out_dir,"figures","Heatmap_Spearman_env_vs_topGenera.png"),
    width = 10, height = 5
  )
}
#
#
# ============================
# GENUS vs contaminants heatmap (Spearman) + PNG export
# ============================

envG <- env_noNi
rownames(envG) <- meta2$SampleID_clean

# Align
genus_corr <- genus_rel[rownames(envG), , drop=FALSE]

# Keep top 30 genera by mean abundance
topG <- 30
top_genera <- names(sort(colMeans(genus_corr), decreasing=TRUE))[1:min(topG, ncol(genus_corr))]
genus_top <- genus_corr[, top_genera, drop=FALSE]

# Spearman correlations: env vars (rows) x genera (cols)
cor_mat_genus <- sapply(colnames(envG), function(v){
  apply(genus_top, 2, function(x){
    suppressWarnings(cor(x, envG[[v]], method="spearman", use="pairwise.complete.obs"))
  })
})
cor_mat_genus <- t(cor_mat_genus)

write.csv(cor_mat_genus, file.path(out_dir,"tables","SpearmanCorr_env_vs_topGenera.csv"))

if(requireNamespace("pheatmap", quietly = TRUE)){
  pheatmap::pheatmap(
    cor_mat_genus,
    cluster_rows=TRUE, cluster_cols=TRUE,
    main=paste0("Spearman: Env vars vs Top ", topG, " Fungal Genera"),
    filename=file.path(out_dir,"figures","Heatmap_Spearman_env_vs_topGenera.png"),
    width=10, height=6
  )
}

#
#
#
#
library(dplyr)
library(tidyr)
library(ggplot2)

# Load significant genera
sig_loc <- read.csv("C:/Users/valer/OneDrive/Desktop/Amplicon_P.R_Sediments/ITS_R/ITS_ampl_sed/IndVal_GENUS_Location_SIG.csv")
sig_genera <- sig_loc$Genus

# Subset genus_rel and convert to data.frame
genus_sig <- genus_rel[, sig_genera, drop=FALSE] %>% 
  as.data.frame()

# Add SampleID + Location (safe alignment)
genus_sig$SampleID <- rownames(genus_sig)

meta_loc <- meta2 %>%
  select(SampleID_clean, Location) %>%
  distinct()

genus_sig <- genus_sig %>%
  left_join(meta_loc, by = c("SampleID" = "SampleID_clean"))

# Pivot to long
genus_long <- genus_sig %>%
  pivot_longer(cols = all_of(sig_genera),
               names_to = "Genus",
               values_to = "RelAbundance")

# Plot
p_box <- ggplot(genus_long,
                aes(x = Location, y = RelAbundance, fill = Location)) +
  geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~Genus, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Significant Fungal Genus Indicators by Location",
       y = "Relative Abundance",
       x = "Location")

# Save
ggsave(file.path(out_dir, "figures", "IndVal_GENUS_Location_Boxplots.png"),
       p_box, width = 12, height = 8, dpi = 300)


# 


envG <- env_noNi
rownames(envG) <- meta2$SampleID_clean

genus_sig_corr <- genus_rel[, sig_genera, drop=FALSE]

cor_mat_sig <- sapply(colnames(envG), function(v){
  apply(genus_sig_corr, 2, function(x){
    suppressWarnings(cor(x, envG[[v]],
                         method="spearman",
                         use="pairwise.complete.obs"))
  })
})
cor_mat_sig <- t(cor_mat_sig)

write.csv(cor_mat_sig,
          file.path(out_dir,"tables",
                    "SpearmanCorr_env_vs_SignificantGenera.csv"))

if(requireNamespace("pheatmap", quietly = TRUE)){
  pheatmap::pheatmap(
    cor_mat_sig,
    cluster_rows=TRUE,
    cluster_cols=TRUE,
    main="Spearman: Metals vs Significant Fungal Genera",
    filename=file.path(out_dir,"figures",
                       "Heatmap_SignificantGenera_vs_Metals.png"),
    width=8,
    height=6
  )
}



table(meta2$Location)
aggregate(env_noNi$`Cu Average`, list(Location=meta2$Location), mean, na.rm=TRUE)
aggregate(env_noNi$`96hrs Survial %`, list(Location=meta2$Location), mean, na.rm=TRUE)



# -----------------------------
# 10) Taxa/ASV vs contaminants: correlation heatmap (ASV-level, filtered)
#     Spearman correlations between ASV rel. abundance and env vars
# -----------------------------
# Use Model A env (no Ni) for n=29 so it includes all sites
envA <- env_noNi
rownames(envA) <- rownames(asv)

# Align rows for correlation
rel_corr <- rel_f[rownames(envA), , drop=FALSE]

# Spearman correlations ASVs x env variables (this can be big; we’ll keep top-variable ASVs)
# Keep top 100 most abundant ASVs for correlation visuals
top_asvs <- names(sort(colSums(rel_corr), decreasing = TRUE))[1:min(100, ncol(rel_corr))]
rel_top <- rel_corr[, top_asvs, drop=FALSE]

# Correlation matrix: ASVs (rows) vs env vars (cols) after transpose
cor_mat <- sapply(colnames(envA), function(v){
  apply(rel_top, 2, function(x) suppressWarnings(cor(x, envA[[v]], method="spearman", use="pairwise.complete.obs")))
})
cor_mat <- t(cor_mat)  # env vars as rows

write.csv(cor_mat, file.path(out_dir,"tables","SpearmanCorr_env_vs_top100ASVs.csv"))

if(has_pheatmap){
  pheatmap::pheatmap(cor_mat,
                     cluster_rows=TRUE, cluster_cols=TRUE,
                     main="Spearman: Env vars vs Top 100 ASVs (rel. abundance)",
                     filename=file.path(out_dir,"figures","Heatmap_Spearman_env_vs_top100ASVs.png"),
                     width=12, height=5)
}

# -----------------------------
message("✅ DONE. Outputs saved in: ", out_dir)



