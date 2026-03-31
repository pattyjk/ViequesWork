# =========================
# 16S from filtered FASTQ -> DADA2 -> Taxonomy -> Phyloseq -> Plots
# =========================

# ---------- Packages ----------
need <- c("dada2","Biostrings","phyloseq","vegan","ggplot2","dplyr","tidyr","readr","stringr","forcats","data.table")
for(p in need){ if(!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org"); library(p, character.only=TRUE) }
theme_set(theme_bw(base_size = 12))

# ---------- Paths (EDIT THESE) ----------
# Use forward slashes on Windows:
filtered_dir <- "C:/Users/valer/OneDrive/Desktop/Amplicon_P.R_Sediments/ViequesWork-main/Sediment/16S_fastq/filtered"
silva_dir    <- "C:/Users/valer/OneDrive/Desktop/Amplicon_P.R_Sediments/SILVA"

# SILVA files (edit the filenames if your copies use different versions)
silva_train   <- file.path(silva_dir, "silva_nr99_v138.1_train_set.fa.gz")
silva_species <- file.path(silva_dir, "silva_species_assignment_v138.1.fa.gz")  # optional

# Metadata (must include Abbreviation and Zone, and a sample key matching FASTQ sample names)
meta_path  <- "Meta data  - Copy of sediment toxicity .csv"   # if this lives elsewhere, set absolute path

# Output
outdir <- "exports"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Options
DO_RAREFY   <- TRUE
RARE_DEPTH  <- NA  # if NA, rarefy to min library size
ABBREV_COL  <- "Abbreviation"   # color by this
ZONE_COL    <- "Zone"           # shape/facet by this
SAMPLE_KEY  <- "SampleID"       # column in metadata used to match sample names

# ---------- Helpers ----------
clean_ids <- function(x){
  x %>%
    str_replace("\\.fastq\\.gz$", "") %>%
    str_replace("\\.fastq$", "") %>%
    str_replace("_R[12]_001$", "") %>%
    str_replace("_R[12]$", "") %>%
    str_replace("_L00[0-9]", "") %>%
    str_replace_all("[^A-Za-z0-9_\\-]", "_") %>%
    str_replace_all("__+", "_") %>%
    str_replace("^_|_$", "")
}

stopifnot(dir.exists(filtered_dir))
stopifnot(file.exists(silva_train))

# ---------- Find FASTQ files; detect paired vs single ----------
fq <- list.files(filtered_dir, pattern = "\\.fastq(\\.gz)?$", full.names = TRUE)
if(length(fq) == 0) stop("No FASTQ files found in: ", filtered_dir)

# Heuristic: paired if both R1 and R2 present
r1 <- fq[grepl("_R1", basename(fq), ignore.case=TRUE)]
r2 <- fq[grepl("_R2", basename(fq), ignore.case=TRUE)]
IS_PAIRED <- length(r1) > 0 && length(r2) > 0

if(IS_PAIRED){
  # harmonize basenames without _R1/_R2
  base_r1 <- clean_ids(basename(r1))
  base_r2 <- clean_ids(basename(r2))
  # Derive sample names by stripping everything after the last "_R"
  sam_r1 <- gsub("(.*)_R1.*","\\1", basename(r1))
  sam_r2 <- gsub("(.*)_R2.*","\\1", basename(r2))
  if(!all(gsub("_R1.*","", basename(r1)) == gsub("_R2.*","", basename(r2)))){
    stop("R1/R2 pairing mismatch — ensure filenames share a common sample prefix before _R1/_R2.")
  }
  samples <- clean_ids(gsub("_R[12].*","", basename(r1)))
  names(r1) <- samples; names(r2) <- samples
} else {
  samples <- clean_ids(gsub("\\.fastq(\\.gz)?$","", basename(fq)))
  names(fq) <- samples
}

# ---------- DADA2 core pipeline ----------
# Quality (optional): plotQualityProfile(r1[1:2]); plotQualityProfile(r2[1:2])
if(IS_PAIRED){
  # Derep
  derepF <- dada2::derepFastq(r1); names(derepF) <- samples
  derepR <- dada2::derepFastq(r2); names(derepR) <- samples
  # Learn errors
  errF <- learnErrors(r1, multithread = TRUE)
  errR <- learnErrors(r2, multithread = TRUE)
  # Denoise
  dadaF <- dada(derepF, err=errF, multithread=TRUE)
  dadaR <- dada(derepR, err=errR, multithread=TRUE)
  # Merge pairs
  mergers <- mergePairs(dadaF, derepF, dadaR, derepR)
  # Sequence table & chimera removal
  seqtab <- makeSequenceTable(mergers)
} else {
  derep <- dada2::derepFastq(fq); names(derep) <- samples
  err   <- learnErrors(fq, multithread = TRUE)
  dada_out <- dada(derep, err=err, multithread=TRUE)
  seqtab <- makeSequenceTable(dada_out)
}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
message("Reads per sample (non-chimeric):")
reads_per_sample <- rowSums(seqtab.nochim)
write.csv(data.frame(Sample=names(reads_per_sample), Reads=reads_per_sample),
          file.path(outdir,"library_sizes_raw.csv"), row.names=FALSE)

# ---------- Assign taxonomy with SILVA ----------
dna <- DNAStringSet(colnames(seqtab.nochim))
names(dna) <- colnames(seqtab.nochim)

# Drop <50 nt (DADA2 cannot classify those)
len <- width(dna)
if(any(len < 50)){
  keep <- which(len >= 50)
  message(sum(len < 50), " ASVs <50 nt removed before taxonomy.")
  dna <- dna[keep]
  seqtab.nochim <- seqtab.nochim[, names(dna), drop=FALSE]
}

taxa <- assignTaxonomy(dna, silva_train, multithread = TRUE)
if(file.exists(silva_species)) {
  taxa <- addSpecies(taxa, silva_species)
}

# ---------- Build phyloseq (add metadata) ----------
otu_mat <- t(seqtab.nochim)  # ASVs as rows
rownames(otu_mat) <- make.unique(rownames(otu_mat))
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

# taxonomy table
canon <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_df <- as.data.frame(taxa)
for(cl in canon) if(!(cl %in% colnames(tax_df))) tax_df[[cl]] <- NA
tax_df <- tax_df[, canon, drop=FALSE]
tax_df <- tax_df[match(rownames(OTU), rownames(tax_df)), , drop=FALSE]
TAX <- tax_table(as.matrix(tax_df))

# metadata 
stopifnot(file.exists(meta_path))
meta <- read.csv(meta_path, check.names = FALSE)
stopifnot(SAMPLE_KEY %in% names(meta), ABBREV_COL %in% names(meta), ZONE_COL %in% names(meta))
# Normalize sample names in metadata to match FASTQ-derived names
meta[[SAMPLE_KEY]] <- clean_ids(meta[[SAMPLE_KEY]])

# Keep only samples present in seqtab
keep_meta <- meta[[SAMPLE_KEY]] %in% colnames(seqtab.nochim)
meta_ps <- meta[keep_meta, , drop=FALSE]
rownames(meta_ps) <- meta_ps[[SAMPLE_KEY]]
# Reorder to match phyloseq sample order
meta_ps <- meta_ps[colnames(seqtab.nochim), , drop=FALSE]
meta_ps[[ABBREV_COL]] <- factor(meta_ps[[ABBREV_COL]])
meta_ps[[ZONE_COL]]   <- factor(meta_ps[[ZONE_COL]])
SAMP <- sample_data(meta_ps)













# ---------- Build phyloseq (add metadata) ----------
# (1) Tax table from DADA2
otu_mat <- t(seqtab.nochim)              # ASVs as rows
rownames(otu_mat) <- make.unique(rownames(otu_mat))
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

canon <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_df <- as.data.frame(taxa)
for(cl in canon) if(!(cl %in% colnames(tax_df))) tax_df[[cl]] <- NA
tax_df <- tax_df[, canon, drop=FALSE]
tax_df <- tax_df[match(rownames(OTU), rownames(tax_df)), , drop=FALSE]
TAX <- tax_table(as.matrix(tax_df))

# (2) Metadata — robust detection & alignment
stopifnot(file.exists(meta_path))
meta <- read.csv(meta_path, check.names = FALSE, na.strings = c("", "NA", "NaN"))

# robust cleaner (no native-pipe gotchas)
clean_ids <- function(x){
  x <- as.character(x)
  x <- gsub("\\.fastq\\.gz$", "", x)
  x <- gsub("\\.fastq$", "", x)
  x <- gsub("_R[12]_001$", "", x)
  x <- gsub("_R[12]$", "", x)
  x <- gsub("_L00[0-9]", "", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
  x <- gsub("__+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

seq_samples_raw   <- rownames(seqtab.nochim)
stopifnot(length(seq_samples_raw) > 0)
seq_samples_clean <- clean_ids(seq_samples_raw)

# score overlap for each metadata column against seq sample names
col_scores <- vapply(names(meta), function(cn){
  vals <- clean_ids(meta[[cn]])
  sum(vals %in% seq_samples_clean)
}, integer(1), USE.NAMES = TRUE)

# prefer Abbreviation if it matches well (>= 80% of samples)
preferred_key <- NA_character_
if("Abbreviation" %in% names(meta)){
  ab_overlap <- sum(clean_ids(meta$Abbreviation) %in% seq_samples_clean)
  if(ab_overlap >= 0.8 * length(seq_samples_clean)) preferred_key <- "Abbreviation"
}

# choose best key
if(!is.na(preferred_key)){
  SAMPLE_KEY <- preferred_key
} else {
  if(length(col_scores) == 0 || max(col_scores) == 0){
    diag <- c(
      "No metadata column matched sequencing samples.",
      sprintf("Detected sequencing samples (first 10): %s", paste(head(seq_samples_clean,10), collapse=", ")),
      "Overlap by metadata column:", paste(sprintf("- %s : %d", names(col_scores), col_scores), collapse="\n")
    )
    writeLines(diag, file.path(outdir, "metadata_matching_report.txt"))
    stop("No matching sample key in metadata. See exports/metadata_matching_report.txt")
  }
  SAMPLE_KEY <- names(col_scores)[which.max(col_scores)]
}
message("Detected sample key column: ", SAMPLE_KEY)

# normalize chosen key and align to sequencing samples
meta[[SAMPLE_KEY]] <- clean_ids(meta[[SAMPLE_KEY]])
idx <- match(seq_samples_clean, meta[[SAMPLE_KEY]])
if(anyNA(idx)){
  missing <- seq_samples_clean[is.na(idx)]
  writeLines(c("Samples missing in metadata:", paste(missing, collapse=", ")),
             file.path(outdir, "metadata_missing_samples.txt"))
}
keep <- !is.na(idx)
meta_ps <- meta[idx[keep], , drop = FALSE]
rownames(meta_ps) <- seq_samples_clean[keep]

# find Abbreviation & Zone even if headers vary
find_col <- function(df, candidates){
  cn <- names(df)
  norm <- tolower(gsub("[^a-z0-9]", "", cn))
  pats <- tolower(gsub("[^a-z0-9]", "", candidates))
  hit <- which(norm %in% pats)
  if(length(hit)) return(cn[hit[1]])
  hit <- which(vapply(pats, function(p) any(grepl(p, norm)), logical(1)))
  if(length(hit)) {
    for(p in pats){
      w <- which(grepl(p, norm))
      if(length(w)) return(cn[w[1]])
    }
  }
  NA_character_
}

ABBREV_COL <- if("Abbreviation" %in% names(meta_ps)) "Abbreviation" else
  find_col(meta_ps, c("Abbreviation","Abbrev","Initials","Individual","Person","ID"))

ZONE_COL <- if("Zone" %in% names(meta_ps)) "Zone" else
  find_col(meta_ps, c("Zone","Region","Area","SiteZone","Stratum","Group","Site"))

if(is.na(ABBREV_COL) || is.na(ZONE_COL)){
  diag <- c(
    sprintf("ABBREV_COL found? %s (%s)", !is.na(ABBREV_COL), ifelse(is.na(ABBREV_COL),"", ABBREV_COL)),
    sprintf("ZONE_COL   found? %s (%s)", !is.na(ZONE_COL),   ifelse(is.na(ZONE_COL),"",   ZONE_COL)),
    "Available columns:", paste(names(meta_ps), collapse=", ")
  )
  writeLines(diag, file.path(outdir, "metadata_missing_columns.txt"))
  stop("Abbreviation or Zone column not found. See exports/metadata_missing_columns.txt")
}

# enforce factors and build sample_data
meta_ps[[ABBREV_COL]] <- factor(meta_ps[[ABBREV_COL]])
meta_ps[[ZONE_COL]]   <- factor(meta_ps[[ZONE_COL]])
SAMP <- sample_data(meta_ps)

# phyloseq object
ps <- phyloseq(OTU, TAX, SAMP)

# expose names for later plotting code
assign("ABBREV_COL", ABBREV_COL, inherits = TRUE)
assign("ZONE_COL",   ZONE_COL,   inherits = TRUE)


# build phyloseq with aligned metadata
ps <- phyloseq(OTU, TAX, SAMP)

# keep these names for the rest of the pipeline
assign("ABBREV_COL", ABBREV_COL, inherits = TRUE)
assign("ZONE_COL",   ZONE_COL,   inherits = TRUE)

ps <- phyloseq(OTU, TAX, SAMP)




# ---- METADATA: robust read + auto-detect sample key / Abbreviation / Zone ----
library(readr)

read_metadata <- function(path){
  stopifnot(file.exists(path))
  # Try CSV first; if it fails, guess delimiter from first line
  df <- tryCatch(readr::read_csv(path, show_col_types = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) {
    first <- readLines(path, n = 1, warn = FALSE)
    delim <- if (grepl("\t", first)) "\t" else if (grepl(";", first)) ";" else ","
    df <- readr::read_delim(path, delim = delim, show_col_types = FALSE)
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  # Drop fully empty rows
  if (nrow(df) > 0) {
    empty_row <- apply(df, 1, function(r) all(is.na(r) | r == ""))
    if (any(empty_row)) df <- df[!empty_row, , drop = FALSE]
  }
  df
}

clean_ids <- function(x){
  x <- as.character(x)
  x <- gsub("\\.fastq\\.gz$", "", x)
  x <- gsub("\\.fastq$", "", x)
  x <- gsub("_R[12]_001$", "", x)
  x <- gsub("_R[12]$", "", x)
  x <- gsub("_L00[0-9]", "", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
  x <- gsub("__+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# 1) Read metadata robustly
meta <- read_metadata(meta_path)
if (nrow(meta) == 0) stop("Metadata file read with 0 rows. Check the delimiter or the path: ", normalizePath(meta_path))

# 2) Sequencing sample names (from your DADA2 table)
seq_samples_raw   <- rownames(seqtab.nochim)
if (length(seq_samples_raw) == 0) stop("No sequencing samples found in seqtab.nochim (rows).")
seq_samples_clean <- clean_ids(seq_samples_raw)

# 3) Pick the best matching metadata column as the sample key
col_scores <- vapply(names(meta), function(cn){
  vals <- clean_ids(meta[[cn]])
  # Guard: if column length mismatches nrow(meta), coerce to length or 0
  if (length(vals) != nrow(meta)) return(0L)
  sum(vals %in% seq_samples_clean)
}, integer(1), USE.NAMES = TRUE)

if (length(col_scores) == 0 || max(col_scores) == 0) {
  writeLines(c(
    "No metadata column matched sequencing samples.",
    sprintf("Sequencing samples (first 10): %s", paste(head(seq_samples_clean, 10), collapse=", ")),
    "Overlap by metadata column:",
    paste(sprintf("- %s : %d", names(col_scores), col_scores), collapse = "\n")
  ), file.path(outdir, "metadata_matching_report.txt"))
  stop("No matching sample key in metadata. See exports/metadata_matching_report.txt")
}

SAMPLE_KEY <- names(col_scores)[which.max(col_scores)]
message("Detected sample key column: ", SAMPLE_KEY)

# 4) Clean & validate the chosen key (must not be zero-length)
tmp <- meta[[SAMPLE_KEY]]
if (length(tmp) == 0) {
  write.csv(utils::head(meta), file.path(outdir, "metadata_preview.csv"), row.names = FALSE)
  stop("Chosen key column has length 0 after read. See exports/metadata_preview.csv and verify your file.")
}
meta[[SAMPLE_KEY]] <- clean_ids(tmp)

# 5) Align metadata rows to the sequencing samples
idx <- match(seq_samples_clean, meta[[SAMPLE_KEY]])
if (anyNA(idx)) {
  missing <- seq_samples_clean[is.na(idx)]
  writeLines(c("Samples missing in metadata:", paste(missing, collapse=", ")),
             file.path(outdir, "metadata_missing_samples.txt"))
}
keep <- !is.na(idx)
meta_ps <- meta[idx[keep], , drop = FALSE]
rownames(meta_ps) <- seq_samples_clean[keep]

# 6) Detect Abbreviation & Zone columns (flexible names)
find_col <- function(df, candidates){
  cn <- names(df)
  norm <- tolower(gsub("[^a-z0-9]", "", cn))
  pats <- tolower(gsub("[^a-z0-9]", "", candidates))
  # exact
  hit <- which(norm %in% pats)
  if (length(hit)) return(cn[hit[1]])
  # contains
  for (p in pats) {
    w <- which(grepl(p, norm))
    if (length(w)) return(cn[w[1]])
  }
  NA_character_
}

ABBREV_COL <- if ("Abbreviation" %in% names(meta_ps)) "Abbreviation" else
  find_col(meta_ps, c("Abbreviation","Abbrev","Initials","Individual","Person","ID"))

ZONE_COL <- if ("Zone" %in% names(meta_ps)) "Zone" else
  find_col(meta_ps, c("Zone","Region","Area","SiteZone","Stratum","Group","Site"))

if (is.na(ABBREV_COL) || is.na(ZONE_COL)) {
  writeLines(c(
    sprintf("ABBREV_COL found? %s (%s)", !is.na(ABBREV_COL), ifelse(is.na(ABBREV_COL),"", ABBREV_COL)),
    sprintf("ZONE_COL   found? %s (%s)", !is.na(ZONE_COL),   ifelse(is.na(ZONE_COL),"",   ZONE_COL)),
    "Available columns:", paste(names(meta_ps), collapse = ", ")
  ), file.path(outdir, "metadata_missing_columns.txt"))
  stop("Abbreviation or Zone column not found. See exports/metadata_missing_columns.txt")
}

# 7) Build sample_data and phyloseq
meta_ps[[ABBREV_COL]] <- factor(meta_ps[[ABBREV_COL]])
meta_ps[[ZONE_COL]]   <- factor(meta_ps[[ZONE_COL]])
SAMP <- sample_data(meta_ps)

ps <- phyloseq(OTU, TAX, SAMP)

# expose so plotting code uses your chosen columns
assign("ABBREV_COL", ABBREV_COL, inherits = TRUE)
assign("ZONE_COL",   ZONE_COL,   inherits = TRUE)






Q   # (exits the Browse> prompt if you're still in it)

# Peek at what we read
dim(meta); names(meta)[1:10]
utils::head(meta[[SAMPLE_KEY]])
readLines(file.path(outdir, "metadata_preview.csv"), n = 5)


# --- Clean sample key & align metadata to sequencing samples ---
# robust cleaner (no base pipe)
clean_ids <- function(x){
  x <- as.character(x)
  x <- gsub("\\.fastq\\.gz$", "", x)
  x <- gsub("\\.fastq$", "", x)
  x <- gsub("_R[12]_001$", "", x)
  x <- gsub("_R[12]$", "", x)
  x <- gsub("_L00[0-9]", "", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
  x <- gsub("__+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# 1) Normalize key in metadata and sequencing sample names
meta[[SAMPLE_KEY]] <- clean_ids(meta[[SAMPLE_KEY]])
seq_samples_clean  <- clean_ids(rownames(seqtab.nochim))
stopifnot(length(seq_samples_clean) > 0)

# 2) Align metadata rows to sequencing samples (same order)
idx <- match(seq_samples_clean, meta[[SAMPLE_KEY]])
if (anyNA(idx)) {
  missing <- seq_samples_clean[is.na(idx)]
  if (length(missing)) {
    writeLines(c("Samples missing in metadata:", paste(missing, collapse=", ")),
               file.path(outdir, "metadata_missing_samples.txt"))
  }
}
keep <- !is.na(idx)
meta_ps <- meta[idx[keep], , drop = FALSE]
rownames(meta_ps) <- seq_samples_clean[keep]

# 3) Lock in Abbreviation (color) and Zone (shape/facet)
#    If Abbreviation is absent, derive it from the sample key (prefix before first "_")
has_abbrev <- "Abbreviation" %in% names(meta_ps)
if (!has_abbrev) {
  meta_ps$Abbreviation <- toupper(sub("^([A-Za-z0-9\\-]+).*", "\\1", meta_ps[[SAMPLE_KEY]]))
  ABBREV_COL <- "Abbreviation"
} else {
  ABBREV_COL <- "Abbreviation"
}

# Try to find Zone; fallback to 'Location'; else guess a categorical column
ZONE_COL <- NA_character_
if ("Zone" %in% names(meta_ps)) {
  ZONE_COL <- "Zone"
} else if ("Location" %in% names(meta_ps)) {
  ZONE_COL <- "Location"
} else {
  # Guess: choose a non-numeric column with 2–8 unique values (not the key or Abbreviation)
  is_cat <- vapply(meta_ps, function(v) is.character(v) || is.factor(v), logical(1))
  cand <- setdiff(names(meta_ps)[is_cat], c(SAMPLE_KEY, ABBREV_COL))
  if (length(cand)) {
    levs <- sapply(cand, function(cn) length(unique(na.omit(meta_ps[[cn]]))))
    cand2 <- cand[levs >= 2 & levs <= 8]
    if (length(cand2)) ZONE_COL <- cand2[1]
  }
}
if (is.na(ZONE_COL)) {
  # last resort: create a placeholder to keep plots running
  meta_ps$Zone <- factor("Unknown")
  ZONE_COL <- "Zone"
  writeLines("ZONE not found; created placeholder 'Zone = Unknown'. Consider mapping from 'Location' or fixing the header.",
             file.path(outdir, "metadata_zone_fallback.txt"))
}

# Ensure factors
meta_ps[[ABBREV_COL]] <- factor(meta_ps[[ABBREV_COL]])
meta_ps[[ZONE_COL]]   <- factor(meta_ps[[ZONE_COL]])

# 4) Build sample_data and (re)build phyloseq
SAMP <- sample_data(meta_ps)
ps   <- phyloseq(OTU, TAX, SAMP)

message("Using columns: Abbreviation = '", ABBREV_COL, "', Zone = '", ZONE_COL, "'.")
# Quick sanity
print(table(meta_ps[[ZONE_COL]]))
print(utils::head(meta_ps[, c(SAMPLE_KEY, ABBREV_COL, ZONE_COL)], 6))

# --- Continue exactly as before (rarefaction, alpha/beta, ordinations) ---

# Rarefaction (optional)
lib_sizes <- sample_sums(ps)
if (DO_RAREFY) {
  if (is.na(RARE_DEPTH)) RARE_DEPTH <- min(lib_sizes)
  set.seed(1)
  ps_use <- suppressWarnings(rarefy_even_depth(ps, sample.size = RARE_DEPTH, rngseed=1, replace=FALSE, verbose=FALSE))
} else {
  ps_use <- ps
}

# Alpha diversity (x = Zone; color = Abbreviation) + export
alpha_df <- estimate_richness(ps_use, measures = c("Shannon","Simpson","Observed"))
alpha_df <- tibble::rownames_to_column(alpha_df, "Sample")
alpha_df <- merge(alpha_df,
                  tibble::rownames_to_column(as.data.frame(meta_ps), "Sample"),
                  by = "Sample", all.x = TRUE)
write.csv(alpha_df, file.path(outdir, "alpha_diversity_byZoneAbbrev.csv"), row.names = FALSE)

for (metric in c("Shannon","Simpson","Observed")) {
  gp <- ggplot(alpha_df, aes(x = .data[[ZONE_COL]], y = .data[[metric]], color = .data[[ABBREV_COL]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.95) +
    labs(title = paste("Alpha diversity:", metric), x = "Zone", y = metric, color = "Individual") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  ggsave(file.path(outdir, paste0("alpha_", tolower(metric), "_byZone_colorAbbrev.png")),
         gp, width = 8, height = 5, dpi = 300)
}

# Beta diversity & ordinations (color = Abbreviation; shape = Zone)
dist_bray <- phyloseq::distance(ps_use, method="bray")
set.seed(1)
ord_nmds <- ordinate(ps_use, method="NMDS", distance=dist_bray, trymax=100)
ord_pcoa <- ordinate(ps_use, method="PCoA", distance=dist_bray)

p_nmds <- plot_ordination(ps_use, ord_nmds, type="samples", color=ABBREV_COL, shape=ZONE_COL) +
  geom_point(size=3, alpha=0.95) +
  labs(title="NMDS (Bray-Curtis): color=Individual, shape=Zone", color="Individual", shape="Zone")
ggsave(file.path(outdir, "ordination_NMDS_bray_colorAbbrev_shapeZone.png"), p_nmds, width=7, height=5, dpi=300)

p_pcoa <- plot_ordination(ps_use, ord_pcoa, type="samples", color=ABBREV_COL, shape=ZONE_COL) +
  geom_point(size=3, alpha=0.95) +
  labs(title="PCoA (Bray-Curtis): color=Individual, shape=Zone", color="Individual", shape="Zone")
ggsave(file.path(outdir, "ordination_PCoA_bray_colorAbbrev_shapeZone.png"), p_pcoa, width=7, height=5, dpi=300)

# Save ordination coordinates
ord_scores <- as.data.frame(scores(ord_nmds, display="sites"))
ord_scores <- tibble::rownames_to_column(ord_scores, "Sample")
ord_scores <- merge(ord_scores,
                    tibble::rownames_to_column(as.data.frame(meta_ps), "Sample"),
                    by="Sample", all.x=TRUE)
write.csv(ord_scores, file.path(outdir, "ordination_NMDS_scores_withZoneAbbrev.csv"), row.names = FALSE)

pcoa_scores <- as.data.frame(ord_pcoa$vectors)
pcoa_scores <- tibble::rownames_to_column(pcoa_scores, "Sample")
pcoa_scores <- merge(pcoa_scores,
                     tibble::rownames_to_column(as.data.frame(meta_ps), "Sample"),
                     by="Sample", all.x=TRUE)
write.csv(pcoa_scores, file.path(outdir, "ordination_PCoA_scores_withZoneAbbrev.csv"), row.names = FALSE)

outdir <- "exports"
path <- normalizePath(outdir, mustWork = TRUE)
# Open File Explorer to that folder:
shell.exec(path)          # Windows
# (Alt) In any OS: utils::browseURL(path)


# List all images we saved
list.files("exports", pattern = "\\.png$", full.names = TRUE)

# Open a specific plot in your default image viewer:
utils::browseURL(normalizePath(file.path("exports", "ordination_NMDS_bray_colorAbbrev_shapeZone.png")))
utils::browseURL(normalizePath(file.path("exports", "alpha_shannon_byZone_colorAbbrev.png")))

## =========================
## PERMANOVA + dispersion (Bray, Jaccard, Unifrac, WUnifrac)
## =========================
suppressPackageStartupMessages({ library(phyloseq); library(vegan); library(dplyr); library(tibble); library(readr) })

stopifnot(inherits(ps_use, "phyloseq"))
md <- as(sample_data(ps_use), "data.frame")
stopifnot(all(c("Zone","Location") %in% names(md)))

# Optional: block permutations if you have repeated measures (e.g., per person)
STRATA_COL <- NULL            # e.g., "Abbreviation" or NULL
N_PERM     <- 999

OUTDIR <- file.path("exports","stats_permanova")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- Distances ----
# Bray–Curtis on relative abundance
ps_rel  <- transform_sample_counts(ps_use, function(x) x/sum(x))
OTU_rel <- as(otu_table(ps_rel), "matrix"); if (taxa_are_rows(ps_rel)) OTU_rel <- t(OTU_rel)
D_bray  <- vegdist(OTU_rel, method = "bray")

# Jaccard (presence/absence)
OTU_pa  <- (as(otu_table(ps_use), "matrix") > 0) * 1; if (taxa_are_rows(ps_use)) OTU_pa <- t(OTU_pa)
D_jac   <- vegdist(OTU_pa, method = "jaccard", binary = TRUE)

# UniFrac needs a tree in ps_tree (same samples)
HAS_TREE <- exists("ps_tree") && inherits(ps_tree, "phyloseq") &&
  !is.null(phy_tree(ps_tree, errorIfNULL = FALSE))
if (!HAS_TREE) message("No phylogenetic tree detected; skipping UniFrac metrics.")

# ---- helpers ----
.adonis_wrap <- function(D, meta, term, strata = NULL, nperm = 999) {
  form <- as.formula(paste("D ~", term))
  adonis2(form, data = meta,
          permutations = nperm,
          strata = if (length(strata) && !is.null(strata)) meta[[strata]] else NULL,
          by = "margin")
}

.save_tbl <- function(x, path) {
  x <- as.data.frame(x); write.csv(x, path, row.names = TRUE)
}

.run_one <- function(D, metric, factor, meta, strata_col = NULL) {
  # PERMANOVA
  ad <- .adonis_wrap(D, meta, term = factor, strata = strata_col, nperm = N_PERM)
  .save_tbl(ad, file.path(OUTDIR, sprintf("permanova_%s_by_%s.csv", metric, factor)))
  
  # Pull the row for the grouping factor
  ad_df <- as.data.frame(ad)
  rname <- if (factor %in% rownames(ad_df)) factor else grep(factor, rownames(ad_df), value = TRUE)[1]
  R2    <- ad_df[rname, "R2"]; Fval <- ad_df[rname, "F"]; pval <- ad_df[rname, "Pr(>F)"]
  
  # Dispersion
  bd <- betadisper(D, meta[[factor]])
  bd_aov  <- as.data.frame(anova(bd))
  bd_perm <- as.data.frame(permutest(bd, permutations = N_PERM)$tab)
  bd_tuk  <- as.data.frame(TukeyHSD(bd)$groups)
  
  .save_tbl(bd_aov,  file.path(OUTDIR, sprintf("betadisper_%s_by_%s_anova.csv",     metric, factor)))
  .save_tbl(bd_perm, file.path(OUTDIR, sprintf("betadisper_%s_by_%s_permutest.csv", metric, factor)))
  .save_tbl(bd_tuk,  file.path(OUTDIR, sprintf("betadisper_%s_by_%s_TukeyHSD.csv",  metric, factor)))
  
  tibble(
    metric = metric, factor = factor,
    permanova_F = as.numeric(Fval),
    permanova_R2 = as.numeric(R2),
    permanova_p = as.numeric(pval),
    dispersion_F = suppressWarnings(as.numeric(bd_perm[1, "F"])),
    dispersion_p = suppressWarnings(as.numeric(bd_perm[1, "Pr(>F)"]))
  )
}

# ---- run everything ----
summary_rows <- list()
for (fac in c("Zone","Location")) {
  summary_rows[[paste0("bray_", fac)]] <- .run_one(D_bray, "BrayCurtis", fac, md, STRATA_COL)
  summary_rows[[paste0("jaccard_", fac)]] <- .run_one(D_jac, "Jaccard", fac, md, STRATA_COL)
  if (HAS_TREE) {
    D_uw <- phyloseq::distance(ps_tree, method = "unifrac")
    D_w  <- phyloseq::distance(ps_tree, method = "wunifrac")
    summary_rows[[paste0("unifrac_unweighted_", fac)]] <- .run_one(D_uw, "UniFrac_unweighted", fac, md, STRATA_COL)
    summary_rows[[paste0("unifrac_weighted_",   fac)]] <- .run_one(D_w,  "UniFrac_weighted",   fac, md, STRATA_COL)
  }
}

summary_df <- dplyr::bind_rows(summary_rows)
write.csv(summary_df, file.path(OUTDIR, "permanova_dispersion_summary.csv"), row.names = FALSE)
message("✓ PERMANOVA + dispersion finished. See: ", normalizePath(OUTDIR))
print(summary_df)


# ==== HEATMAPS ====
# ==== HEATMAPS (robust numeric coercion + safe alignment) ====
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap", repos="https://cloud.r-project.org")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", repos="https://cloud.r-project.org")
library(pheatmap); library(RColorBrewer)

# Helper: force numeric matrix (avoids 'cut.default ... x must be numeric')
num_mat <- function(m){
  m <- as.matrix(m)
  suppressWarnings(storage.mode(m) <- "numeric")   # force numeric
  m[is.na(m)] <- 0
  m
}

# Use the detected column names (set earlier by the script)
abbrev_col <- if (exists("ABBREV_COL")) ABBREV_COL else "Abbreviation"
zone_col   <- if (exists("ZONE_COL"))   ZONE_COL   else "Zone"

# Sample annotations
sdat <- as.data.frame(sample_data(ps_use))
stopifnot(all(c(zone_col, abbrev_col) %in% names(sdat)))
ann  <- sdat[, c(zone_col, abbrev_col), drop = FALSE]
ann[[zone_col]]   <- factor(ann[[zone_col]])
ann[[abbrev_col]] <- factor(ann[[abbrev_col]])
rownames(ann)     <- rownames(sdat)

# Order columns by Zone then Abbreviation
ord_cols <- order(ann[[zone_col]], ann[[abbrev_col]])
ann      <- ann[ord_cols, , drop = FALSE]

# Annotation colors (categorical)
zones <- levels(ann[[zone_col]])
abbrs <- levels(ann[[abbrev_col]])

zone_cols <- setNames(brewer.pal(max(3, min(8, length(zones))), "Set2")[seq_along(zones)], zones)
make_big_palette <- function(n){ hues <- seq(15, 375, length.out = n + 1); hcl(h=hues, l=65, c=100)[1:n] }
abbr_cols <- setNames(make_big_palette(length(abbrs)), abbrs)
ann_colors <- list()
ann_colors[[zone_col]]   <- zone_cols
ann_colors[[abbrev_col]] <- abbr_cols

# Safe aligner for heatmap matrices
align_to_ann <- function(mat){
  mat <- num_mat(mat)
  common <- intersect(colnames(mat), rownames(ann))
  if (length(common) == 0) stop("No overlapping sample names between matrix and annotations.")
  mat <- mat[, common, drop = FALSE]
  ann <<- ann[common, , drop = FALSE]   # update global 'ann' to the same order
  mat
}

# Prefer Genus if taxonomy present; else ASV
has_genus <- !is.null(tax_table(ps_use)) &&
  "Genus" %in% colnames(tax_table(ps_use)) &&
  !all(is.na(tax_table(ps_use)[, "Genus"]))

if (has_genus) {
  # ---- Genus heatmaps ----
  ps_genus   <- tax_glom(ps_use, taxrank = "Genus", NArm = TRUE)
  ps_genus_r <- transform_sample_counts(ps_genus, function(x) x/sum(x))  # relative abundance
  mat <- as.matrix(otu_table(ps_genus_r))
  rownames(mat) <- as.vector(tax_table(ps_genus_r)[, "Genus"])
  rownames(mat)[is.na(rownames(mat))] <- "Unclassified"
  # Collapse duplicate genus names
  mat <- rowsum(mat, group = rownames(mat))
  
  # Top 30 genera
  topN <- 30
  top_rows <- order(rowMeans(mat), decreasing = TRUE)[1:min(topN, nrow(mat))]
  mat_top <- mat[top_rows, , drop = FALSE]
  mat_top <- align_to_ann(mat_top)
  
  # 1) Row z-score
  mat_z <- t(scale(t(mat_top))); mat_z <- num_mat(mat_z)
  zcols <- colorRampPalette(c("#073b4c", "#ffd166", "#ef476f"))(101)
  png(file.path(outdir, "heatmap_genus_top30_zscore.png"), width = 1600, height = 900, res = 140)
  pheatmap(
    mat_z,
    annotation_col    = ann,
    annotation_colors = ann_colors,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    cluster_cols      = FALSE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main              = "Top 30 genera • row z-score of relative abundance",
    color             = zcols
  )
  dev.off()
  
  # 2) Log-scaled relative abundance
  mat_log <- log10(mat_top * 1000 + 1); mat_log <- num_mat(mat_log)
  png(file.path(outdir, "heatmap_genus_top30_log.png"), width = 1600, height = 900, res = 140)
  pheatmap(
    mat_log,
    annotation_col    = ann,
    annotation_colors = ann_colors,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    cluster_cols      = FALSE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main              = "Top 30 genera • log10(relative abundance × 1000 + 1)"
  )
  dev.off()
  
  write.csv(mat_top, file.path(outdir, "heatmap_genus_top30_relabund_matrix.csv"))
  write.csv(mat_z,   file.path(outdir, "heatmap_genus_top30_rowscaled_matrix.csv"))
  
} else {
  # ---- ASV heatmaps (no taxonomy) ----
  ps_rel <- transform_sample_counts(ps_use, function(x) x/sum(x))
  mat <- as.matrix(otu_table(ps_rel))                       # rows = ASVs, cols = samples
  topN <- 50
  top_rows <- order(rowMeans(mat), decreasing = TRUE)[1:min(topN, nrow(mat))]
  mat_top <- mat[top_rows, , drop = FALSE]
  mat_top <- align_to_ann(mat_top)
  
  # Row z-score
  mat_z <- t(scale(t(mat_top))); mat_z <- num_mat(mat_z)
  zcols <- colorRampPalette(c("#073b4c", "#ffd166", "#ef476f"))(101)
  png(file.path(outdir, "heatmap_ASV_top50_zscore.png"), width = 1600, height = 900, res = 140)
  pheatmap(
    mat_z,
    annotation_col    = ann,
    annotation_colors = ann_colors,
    show_rownames     = FALSE,   # ASV IDs are long
    show_colnames     = FALSE,
    cluster_cols      = FALSE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main              = "Top 50 ASVs • row z-score of relative abundance",
    color             = zcols
  )
  dev.off()
  
  # Log-scaled
  mat_log <- log10(mat_top * 1000 + 1); mat_log <- num_mat(mat_log)
  png(file.path(outdir, "heatmap_ASV_top50_log.png"), width = 1600, height = 900, res = 140)
  pheatmap(
    mat_log,
    annotation_col    = ann,
    annotation_colors = ann_colors,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    cluster_cols      = FALSE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main              = "Top 50 ASVs • log10(relative abundance × 1000 + 1)"
  )
  dev.off()
  
  write.csv(mat_top, file.path(outdir, "heatmap_ASV_top50_relabund_matrix.csv"))
  write.csv(mat_z,   file.path(outdir, "heatmap_ASV_top50_rowscaled_matrix.csv"))
}




utils::browseURL(normalizePath(file.path("exports","heatmap_genus_top30_zscore.png")))











# ==== HEATMAPS (robust numeric coercion + safe alignment) ====
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap", repos="https://cloud.r-project.org")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", repos="https://cloud.r-project.org")
library(pheatmap); library(RColorBrewer)

# Helper: force numeric matrix (avoids 'cut.default ... x must be numeric')
num_mat <- function(m){
  m <- as.matrix(m)
  suppressWarnings(storage.mode(m) <- "numeric")   # force numeric
  m[is.na(m)] <- 0
  m
}

# Use the detected column names (set earlier by the script)
abbrev_col <- if (exists("ABBREV_COL")) ABBREV_COL else "Abbreviation"
zone_col   <- if (exists("ZONE_COL"))   ZONE_COL   else "Zone"

# Sample annotations
sdat <- as.data.frame(sample_data(ps_use))
stopifnot(all(c(zone_col, abbrev_col) %in% names(sdat)))
ann  <- sdat[, c(zone_col, abbrev_col), drop = FALSE]
ann[[zone_col]]   <- factor(ann[[zone_col]])
ann[[abbrev_col]] <- factor(ann[[abbrev_col]])
rownames(ann)     <- rownames(sdat)

# Order columns by Zone then Abbreviation
ord_cols <- order(ann[[zone_col]], ann[[abbrev_col]])
ann      <- ann[ord_cols, , drop = FALSE]

# Annotation colors (categorical)
zones <- levels(ann[[zone_col]])
abbrs <- levels(ann[[abbrev_col]])

zone_cols <- setNames(brewer.pal(max(3, min(8, length(zones))), "Set2")[seq_along(zones)], zones)
make_big_palette <- function(n){ hues <- seq(15, 375, length.out = n + 1); hcl(h=hues, l=65, c=100)[1:n] }
abbr_cols <- setNames(make_big_palette(length(abbrs)), abbrs)
ann_colors <- list()
ann_colors[[zone_col]]   <- zone_cols
ann_colors[[abbrev_col]] <- abbr_cols

# Safe aligner for heatmap matrices
align_to_ann <- function(mat){
  mat <- num_mat(mat)
  common <- intersect(colnames(mat), rownames(ann))
  if (length(common) == 0) stop("No overlapping sample names between matrix and annotations.")
  mat <- mat[, common, drop = FALSE]
  ann <<- ann[common, , drop = FALSE]   # update global 'ann' to the same order
  mat
}

# Prefer Genus if taxonomy present; else ASV
has_genus <- !is.null(tax_table(ps_use)) &&
             "Genus" %in% colnames(tax_table(ps_use)) &&
             !all(is.na(tax_table(ps_use)[, "Genus"]))

if (has_genus) {
  # ---- Genus heatmaps ----
  ps_genus   <- tax_glom(ps_use, taxrank = "Genus", NArm = TRUE)
  ps_genus_r <- transform_sample_counts(ps_genus, function(x) x/sum(x))  # relative abundance
  mat <- as.matrix(otu_table(ps_genus_r))
  rownames(mat) <- as.vector(tax_table(ps_genus_r)[, "Genus"])
  rownames(mat)[is.na(rownames(mat))] <- "Unclassified"
  # Collapse duplicate genus names
  mat <- rowsum(mat, group = rownames(mat))

  # Top 30 genera
  topN <- 30
  top_rows <- order(rowMeans(mat), decreasing = TRUE)[1:min(topN, nrow(mat))]
  mat_top <- mat[top_rows, , drop = FALSE]
  mat_top <- align_to_ann(mat_top)

  # 1) Row z-score
  mat_z <- t(scale(t(mat_top))); mat_z <- num_mat(mat_z)
  zcols <- colorRampPalette(c("#073b4c", "#ffd166", "#ef476f"))(101)
  png(file.path(outdir, "heatmap_genus_top30_zscore.png"), width = 1600, height = 900, res = 140)
  pheatmap(
    mat_z,
    annotation_col    = ann,
    annotation_colors = ann_colors,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    cluster_cols      = FALSE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main              = "Top 30 genera • row z-score of relative abundance",
    color             = zcols
  )
  dev.off()

  # 2) Log-scaled relative abundance
  mat_log <- log10(mat_top * 1000 + 1); mat_log <- num_mat(mat_log)
  png(file.path(outdir, "heatmap_genus_top30_log.png"), width = 1600, height = 900, res = 140)
  pheatmap(
    mat_log,
    annotation_col    = ann,
    annotation_colors = ann_colors,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    cluster_cols      = FALSE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main              = "Top 30 genera • log10(relative abundance × 1000 + 1)"
  )
  dev.off()

  write.csv(mat_top, file.path(outdir, "heatmap_genus_top30_relabund_matrix.csv"))
  write.csv(mat_z,   file.path(outdir, "heatmap_genus_top30_rowscaled_matrix.csv"))

} else {
  # ---- ASV heatmaps (no taxonomy) ----
  ps_rel <- transform_sample_counts(ps_use, function(x) x/sum(x))
  mat <- as.matrix(otu_table(ps_rel))                       # rows = ASVs, cols = samples
  topN <- 50
  top_rows <- order(rowMeans(mat), decreasing = TRUE)[1:min(topN, nrow(mat))]
  mat_top <- mat[top_rows, , drop = FALSE]
  mat_top <- align_to_ann(mat_top)

  # Row z-score
  mat_z <- t(scale(t(mat_top))); mat_z <- num_mat(mat_z)
  zcols <- colorRampPalette(c("#073b4c", "#ffd166", "#ef476f"))(101)
  png(file.path(outdir, "heatmap_ASV_top50_zscore.png"), width = 1600, height = 900, res = 140)
  pheatmap(
    mat_z,
    annotation_col    = ann,
    annotation_colors = ann_colors,
    show_rownames     = FALSE,   # ASV IDs are long
    show_colnames     = FALSE,
    cluster_cols      = FALSE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main              = "Top 50 ASVs • row z-score of relative abundance",
    color             = zcols
  )
  dev.off()

  # Log-scaled
  mat_log <- log10(mat_top * 1000 + 1); mat_log <- num_mat(mat_log)
  png(file.path(outdir, "heatmap_ASV_top50_log.png"), width = 1600, height = 900, res = 140)
  pheatmap(
    mat_log,
    annotation_col    = ann,
    annotation_colors = ann_colors,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    cluster_cols      = FALSE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main              = "Top 50 ASVs • log10(relative abundance × 1000 + 1)"
  )
  dev.off()

  write.csv(mat_top, file.path(outdir, "heatmap_ASV_top50_relabund_matrix.csv"))
  write.csv(mat_z,   file.path(outdir, "heatmap_ASV_top50_rowscaled_matrix.csv"))
}

utils::browseURL(normalizePath(file.path("exports","heatmap_genus_top30_zscore.png")))



# If you already built ps/ps_use earlier:
has_tax <- !is.null(tax_table(ps)) && nrow(tax_table(ps)) > 0
all_na  <- if (has_tax) all(is.na(tax_table(ps)[,1])) else TRUE

cat("Has TAX table? ", has_tax, "\nAll NAs in Kingdom? ", all_na, "\n")
if (has_tax && !all_na) head(as.data.frame(tax_table(ps)))



# --- Edit these to your actual files ---
silva_train   <- "C:/Users/valer/OneDrive/Desktop/Amplicon_P.R_Sediments/SILVA/silva_nr99_v138.1_train_set.fa.gz"
silva_species <- "C:/Users/valer/OneDrive/Desktop/Amplicon_P.R_Sediments/SILVA/silva_species_assignment_v138.1.fa.gz"  # optional

stopifnot(file.exists(silva_train))
HAVE_SPECIES <- file.exists(silva_species)

# Get ASV sequences from your table
library(Biostrings); library(dada2); library(phyloseq)
asv_seqs <- colnames(seqtab.nochim)
dna <- DNAStringSet(asv_seqs); names(dna) <- asv_seqs

# Drop <50nt (DADA2 can’t classify those)
len <- width(dna)
if (any(len < 50)) {
  dna <- dna[len >= 50]
  seqtab.nochim <- seqtab.nochim[, names(dna), drop = FALSE]
  message(sum(len < 50), " ASVs <50 nt removed prior to taxonomy.")
}

# Assign taxonomy
tax <- assignTaxonomy(dna, silva_train, multithread = TRUE, minBoot = 50)
if (HAVE_SPECIES) tax <- addSpecies(tax, silva_species)

# Build TAX aligned to your OTU/ASV order
otu_mat <- t(seqtab.nochim)                         # ASVs as rows
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

canon <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_df <- as.data.frame(tax)
for (cl in canon) if (!(cl %in% names(tax_df))) tax_df[[cl]] <- NA
tax_df <- tax_df[match(rownames(OTU), rownames(tax_df)), canon, drop = FALSE]

TAX <- tax_table(as.matrix(tax_df))

# If you already had sample_data, reattach it; else create minimal sample_data
if (exists("meta_ps")) {
  SAMP <- sample_data(meta_ps)
} else {
  SAMP <- sample_data(data.frame(row.names = colnames(seqtab.nochim)))
}

ps <- phyloseq(OTU, TAX, SAMP)

# Save taxonomy for your records
if (!dir.exists("exports")) dir.create("exports")
write.csv(tax_df, "exports/assigned_taxonomy_SILVA.csv")


tx <- as.data.frame(tax_table(ps))
# % classified at each rank
rate <- sapply(tx, function(col) mean(!is.na(col)))
print(round(rate, 3))  # e.g., Kingdom ~1.00, Phylum ~0.98, Genus ~0.7 etc.

# Count of "unclassified" (NA) by rank
na_counts <- sapply(tx, function(col) sum(is.na(col)))
print(na_counts)

# Peek at the most abundant taxa to ensure they look reasonable
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))
top_phylum <- sort(taxa_sums(tax_glom(ps_rel, "Phylum")), decreasing = TRUE)[1:10]
print(top_phylum)

# Sanity: remove non-target sequences (optional)
is_bad <- !is.na(tx$Order) & tx$Order %in% c("Chloroplast","Mitochondria")
sum(is_bad)  # how many chloroplast/mito hits?
if (sum(is_bad) > 0) {
  ps <- prune_taxa(!is_bad, ps)
  message("Removed chloroplast/mitochondria ASVs: ", sum(is_bad))
}




if (!requireNamespace("DECIPHER", quietly = TRUE)) install.packages("DECIPHER", repos="https://cloud.r-project.org")
library(DECIPHER)

# Download/point to DECIPHER’s SILVA training set if you have it (one-time); example path:
# trainingSet <- "C:/.../SILVA_SSU_r138_2021.RData"
# load(trainingSet)  # loads object "trainingSet"

# If you already have DNAStringSet 'dna' aligned to your current ASVs:
ids <- IdTaxa(dna, trainingSet, strand="both", processors=1)
# Convert to taxonomy table
id_to_df <- function(ids){
  ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus")
  out <- matrix(NA_character_, nrow=length(ids), ncol=length(ranks), dimnames=list(names(ids), ranks))
  for (i in seq_along(ids)){
    lab <- ids[[i]]$taxon
    rk  <- sub(":.*$", "", lab)
    nm  <- sub("^[^:]+: ?", "", lab)
    keep <- rk %in% ranks
    out[i, rk[keep]] <- nm[keep]
  }
  as.data.frame(out)
}
tax_df_id <- id_to_df(ids)

# Compare to DADA2 at Genus
common <- intersect(rownames(tax_df_id), rownames(tax_df))
table(DADA2 = tax_df[common,"Genus"], IdTaxa = tax_df_id[common,"Genus"])[1:10, 1:10]






# Phylum barplot (top 12)
ps_phylum <- tax_glom(ps, "Phylum", NArm=TRUE)
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x/sum(x))
top12 <- names(sort(taxa_sums(ps_phylum_rel), decreasing=TRUE))[1:min(12, ntaxa(ps_phylum_rel))]
df <- psmelt(prune_taxa(top12, ps_phylum_rel))

library(ggplot2)
ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_col() + labs(title="Relative abundance by Phylum", y="Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



##taxanomy attached 
stopifnot(inherits(ps, "phyloseq"))
tx  <- tax_table(ps)
otu <- otu_table(ps)

cat("Tax rows:", nrow(tx), "  OTU taxa:", ntaxa(ps), "\n")
cat("Row alignment (taxonomy rows == taxa_names(ps))? ",
    identical(rownames(tx), taxa_names(ps)), "\n")

# Quick look at ranks
head(as.data.frame(tx)[, c("Kingdom","Phylum","Class","Order","Family","Genus","Species")])




txdf <- as.data.frame(tx)
coverage <- sapply(txdf[c("Kingdom","Phylum","Class","Order","Family","Genus","Species")],
                   function(col) mean(!is.na(col)))
round(coverage, 3)

txdf <- as.data.frame(tax_table(ps))
is_euk  <- !is.na(txdf$Kingdom) & txdf$Kingdom != "Bacteria"
is_chl  <- !is.na(txdf$Order)   & txdf$Order == "Chloroplast"
is_mito <- !is.na(txdf$Family)  & txdf$Family == "Mitochondria"
bad <- is_euk | is_chl | is_mito

message("Removing ", sum(bad, na.rm=TRUE), " non-target ASVs (Euk/Chloroplast/Mitochondria).")
ps_clean <- prune_taxa(!bad, ps)




old_ids <- taxa_names(ps_clean)               # sequences
new_ids <- paste0("ASV", seq_along(old_ids))  # ASV1..N

# add the sequence as a new taxonomy column
TAX <- tax_table(ps_clean)
TAX2 <- cbind(as.matrix(TAX), Sequence = old_ids)
tax_table(ps_clean) <- tax_table(TAX2)

# now rename taxa to short IDs
taxa_names(ps_clean) <- new_ids

# export a map for your records
write.csv(data.frame(ASV=new_ids, Sequence=old_ids,
                     as.data.frame(tax_table(ps_clean))),
          file = file.path("exports","asv_id_map_with_taxonomy.csv"),
          row.names = FALSE)


# Phylum composition (relative)
ps_phylum <- tax_glom(ps_clean, "Phylum", NArm=TRUE)
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x/sum(x))
sort(taxa_sums(ps_phylum_rel), decreasing=TRUE)[1:10]

# Genus heatmap (the block I gave before will now run the *genus* path automatically)

# Force any object -> numeric matrix (keeps dimnames)
to_num_mat <- function(m){
  m <- as.matrix(m)
  suppressWarnings(storage.mode(m) <- "double")
  m[!is.finite(m)] <- 0
  m[is.na(m)] <- 0
  m
}

# Align a matrix's columns to your sample annotations
align_to_ann <- function(mat, ann){
  mat <- to_num_mat(mat)
  common <- intersect(colnames(mat), rownames(ann))
  if (length(common) == 0) stop("No overlapping sample names between heatmap matrix and annotations.")
  mat <- mat[, common, drop = FALSE]
  ann <- ann[common, , drop = FALSE]
  list(mat = mat, ann = ann)
}


# Use the detected names set earlier in your script
abbrev_col <- if (exists("ABBREV_COL")) ABBREV_COL else "Abbreviation"
zone_col   <- if (exists("ZONE_COL"))   ZONE_COL   else "Zone"

sdat <- as.data.frame(sample_data(ps_use))
stopifnot(all(c(zone_col, abbrev_col) %in% names(sdat)))
ann  <- sdat[, c(zone_col, abbrev_col), drop = FALSE]
ann[[zone_col]]   <- factor(ann[[zone_col]])
ann[[abbrev_col]] <- factor(ann[[abbrev_col]])
rownames(ann)     <- rownames(sdat)

# Sort columns by Zone then Abbreviation for readability
ord <- order(ann[[zone_col]], ann[[abbrev_col]])
ann <- ann[ord, , drop = FALSE]

# Simple discrete palettes
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", repos="https://cloud.r-project.org")
library(RColorBrewer)
zones <- levels(ann[[zone_col]])
abbrs <- levels(ann[[abbrev_col]])
zone_cols <- setNames(brewer.pal(max(3, min(8, length(zones))), "Set2")[seq_along(zones)], zones)
make_big_palette <- function(n){ hues <- seq(15, 375, length.out = n + 1); hcl(h=hues, l=65, c=100)[1:n] }
abbr_cols <- setNames(make_big_palette(length(abbrs)), abbrs)
ann_colors <- list(); ann_colors[[zone_col]] <- zone_cols; ann_colors[[abbrev_col]] <- abbr_cols




if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap", repos="https://cloud.r-project.org")
library(pheatmap)

has_genus <- !is.null(tax_table(ps_use)) &&
  "Genus" %in% colnames(tax_table(ps_use)) &&
  !all(is.na(tax_table(ps_use)[, "Genus"]))

if (has_genus) {
  # Genus matrix (relative abundance)
  ps_genus   <- tax_glom(ps_use, taxrank = "Genus", NArm = TRUE)
  ps_genus_r <- transform_sample_counts(ps_genus, function(x) x / sum(x))
  mat <- as.matrix(otu_table(ps_genus_r))
  rownames(mat) <- as.vector(tax_table(ps_genus_r)[, "Genus"])
  rownames(mat)[is.na(rownames(mat))] <- "Unclassified"
  mat <- rowsum(mat, group = rownames(mat))        # collapse duplicate genus names
  
  # Top genera by mean abundance
  topN <- 30
  top_rows <- order(rowMeans(mat), decreasing = TRUE)[1:min(topN, nrow(mat))]
  mat_top <- mat[top_rows, , drop = FALSE]
  aligned <- align_to_ann(mat_top, ann); mat_top <- aligned$mat; ann <- aligned$ann
  
  # 1) Row z-score heatmap
  mat_z <- t(scale(t(mat_top))); mat_z <- to_num_mat(mat_z)
  cols_z <- colorRampPalette(c("#073b4c", "#ffd166", "#ef476f"))(101)
  png(file.path(outdir, "heatmap_genus_top30_zscore.png"), width = 1600, height = 900, res = 140)
  pheatmap(
    mat_z,
    annotation_col = ann, annotation_colors = ann_colors,
    show_rownames = TRUE, show_colnames = FALSE,
    cluster_cols = FALSE, clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main = "Top 30 genera • row z-score of relative abundance",
    color = cols_z
  ); dev.off()
  
  # 2) Log-scaled abundance heatmap
  mat_log <- log10(mat_top * 1000 + 1); mat_log <- to_num_mat(mat_log)
  png(file.path(outdir, "heatmap_genus_top30_log.png"), width = 1600, height = 900, res = 140)
  pheatmap(
    mat_log,
    annotation_col = ann, annotation_colors = ann_colors,
    show_rownames = TRUE, show_colnames = FALSE,
    cluster_cols = FALSE, clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main = "Top 30 genera • log10(relative abundance × 1000 + 1)"
  ); dev.off()
  
  # Export matrices
  write.csv(mat_top, file.path(outdir, "heatmap_genus_top30_relabund_matrix.csv"))
  write.csv(mat_z,   file.path(outdir, "heatmap_genus_top30_rowscaled_matrix.csv"))
} else {
  message("No usable Genus taxonomy detected — skipping genus heatmap. (You can re-run after taxonomy assignment.)")
}


stopifnot(is.numeric(mat_z[1,1]))
stopifnot(identical(colnames(mat_z), rownames(ann)))


utils::browseURL(normalizePath(file.path("exports","heatmap_genus_top30_zscore.png")))
utils::browseURL(normalizePath(file.path("exports","heatmap_genus_top30_log.png")))
















utils::browseURL(normalizePath(file.path("exports","heatmap_genus_top30_zscore.png")))



































# ---------- Rarefaction (optional) ----------
lib_sizes <- sample_sums(ps)
write.csv(data.frame(Sample=names(lib_sizes), Reads=as.integer(lib_sizes)),
          file.path(outdir, "library_sizes_after_denoise.csv"), row.names = FALSE)
if(DO_RAREFY){
  if(is.na(RARE_DEPTH)) RARE_DEPTH <- min(lib_sizes)
  set.seed(1)
  ps_use <- suppressWarnings(rarefy_even_depth(ps, sample.size = RARE_DEPTH, rngseed=1, replace=FALSE, verbose=FALSE))
} else {
  ps_use <- ps
}


# ================== GENUS HEATMAP (robust numeric) ==================
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap", repos="https://cloud.r-project.org")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", repos="https://cloud.r-project.org")
library(pheatmap)
library(RColorBrewer)

# Annotations (color = Abbreviation; group = Zone) ----
abbrev_col <- if (exists("ABBREV_COL")) ABBREV_COL else "Abbreviation"
zone_col   <- if (exists("ZONE_COL"))   ZONE_COL   else "Zone"

sdat <- as.data.frame(sample_data(ps_use))
stopifnot(all(c(zone_col, abbrev_col) %in% names(sdat)))
ann <- sdat[, c(zone_col, abbrev_col), drop = FALSE]
ann[[zone_col]]   <- factor(ann[[zone_col]])
ann[[abbrev_col]] <- factor(ann[[abbrev_col]])
rownames(ann)     <- rownames(sdat)

# order columns by Zone then Abbreviation (just for readability)
ord <- order(ann[[zone_col]], ann[[abbrev_col]])
ann <- ann[ord, , drop = FALSE]

# discrete palettes
zones <- levels(ann[[zone_col]])
abbrs <- levels(ann[[abbrev_col]])
zone_cols <- setNames(brewer.pal(max(3, min(8, length(zones))), "Set2")[seq_along(zones)], zones)
make_big_palette <- function(n){ hues <- seq(15, 375, length.out = n + 1); hcl(h=hues, l=65, c=100)[1:n] }
abbr_cols <- setNames(make_big_palette(length(abbrs)), abbrs)
ann_colors <- list(); ann_colors[[zone_col]] <- zone_cols; ann_colors[[abbrev_col]] <- abbr_cols

# Build genus matrix (relative abundance) ----
ps_genus   <- tax_glom(ps_use, taxrank = "Genus", NArm = TRUE)
ps_genus_r <- transform_sample_counts(ps_genus, function(x) x / sum(x))

mat <- as.matrix(otu_table(ps_genus_r))
if (!taxa_are_rows(ps_genus_r)) mat <- t(mat)    # ensure rows = taxa, cols = samples
storage.mode(mat) <- "double"                    # FORCE NUMERIC

# Use genus names as rownames; collapse duplicates
g_names <- as.vector(tax_table(ps_genus_r)[, "Genus"])
g_names[is.na(g_names) | g_names == ""] <- "Unclassified"
rownames(mat) <- g_names
mat <- rowsum(mat, group = rownames(mat))        # collapse by identical genus name

# choose top genera by mean rel. abundance
topN <- 30
top_rows <- head(order(rowMeans(mat), decreasing = TRUE), n = min(topN, nrow(mat)))
mat_top <- mat[top_rows, , drop = FALSE]

# align columns to annotation samples (and same order)
common <- intersect(colnames(mat_top), rownames(ann))
if (length(common) == 0) stop("No overlapping sample names between heatmap matrix and annotations.")
mat_top <- mat_top[, common, drop = FALSE]
ann     <- ann[common, , drop = FALSE]

# z-score by row; replace NA/Inf; ensure numeric
mat_z <- t(scale(t(mat_top)))
mat_z[!is.finite(mat_z)] <- 0
storage.mode(mat_z) <- "double"                  # FORCE NUMERIC

# sanity checks (should pass)
stopifnot(is.numeric(mat_top[1,1]), is.numeric(mat_z[1,1]))
stopifnot(identical(colnames(mat_z), rownames(ann)))

# color palette for z-score heatmap
cols_z <- colorRampPalette(c("#073b4c", "#ffd166", "#ef476f"))(101)

# write z-score heatmap
png(file.path(outdir, "heatmap_genus_top30_zscore.png"), width = 1600, height = 900, res = 140)
pheatmap(
  mat_z,
  annotation_col    = ann,
  annotation_colors = ann_colors,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  cluster_cols      = FALSE,     # keep our Zone->Abbrev order
  clustering_method = "complete",
  clustering_distance_rows = "euclidean",
  main              = "Top 30 genera • row z-score of relative abundance",
  color             = cols_z
)
dev.off()

# Optional: log-scaled heatmap (comment this block out if you don't need it)
mat_log <- log10(mat_top * 1000 + 1)
mat_log[!is.finite(mat_log)] <- 0
storage.mode(mat_log) <- "double"


str(mat_z[1:3, 1:3])
sapply(mat_z[,1:min(5,ncol(mat_z)), drop=FALSE], class)




png(file.path(outdir, "heatmap_genus_top30_log.png"), width = 1600, height = 900, res = 140)
pheatmap(
  mat_log,
  annotation_col    = ann,
  annotation_colors = ann_colors,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  cluster_cols      = FALSE,
  clustering_method = "complete",
  clustering_distance_rows = "euclidean",
  main              = "Top 30 genera • log10(relative abundance × 1000 + 1)"
)
dev.off()

# export matrices for Excel/archival
write.csv(mat_top, file.path(outdir, "heatmap_genus_top30_relabund_matrix.csv"))
write.csv(mat_z,   file.path(outdir, "heatmap_genus_top30_rowscaled_matrix.csv"))
# ================================================================


str(mat_z[1:3, 1:3])
str(ann[1:3, , drop = FALSE])
dput(head(colnames(mat_z)))
dput(head(rownames(ann)))




# ---- make sure 'mat_z' and 'ann' exist as in your current code ----
stopifnot(exists("mat_z"), is.matrix(mat_z))
stopifnot(exists("ann"), is.data.frame(ann))
# force numeric one more time and align
storage.mode(mat_z) <- "double"
mat_z[!is.finite(mat_z)] <- 0
mat_z <- mat_z[, rownames(ann), drop = FALSE]

# explicit breaks for z-scores
ncols <- 101
cols_z <- colorRampPalette(c("#073b4c", "#ffd166", "#ef476f"))(ncols)
rng <- range(mat_z, na.rm = TRUE)
if (diff(rng) < .Machine$double.eps) rng <- rng + c(-1e-6, 1e-6)  # avoid zero-range
brks <- seq(rng[1], rng[2], length.out = ncols + 1)

# try pheatmap, fall back to ggplot if it errors
ok <- TRUE
tryCatch({
  png(file.path(outdir, "heatmap_genus_top30_zscore.png"), width = 1600, height = 900, res = 140)
  pheatmap::pheatmap(
    mat_z,
    color = cols_z,
    breaks = brks,
    annotation_col    = ann,
    annotation_colors = ann_colors,  # from previous block
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    cluster_cols      = FALSE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    main              = "Top 30 genera • row z-score of relative abundance"
  )
  dev.off()
}, error = function(e){ ok <<- FALSE; message("pheatmap failed: ", conditionMessage(e)) })




if (!ok) {
  message("Falling back to ggplot2 heatmap…")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos="https://cloud.r-project.org")
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr", repos="https://cloud.r-project.org")
  library(ggplot2); library(tidyr)
  
  # long format
  df_long <- as.data.frame(mat_z)
  df_long$Genus <- rownames(mat_z)
  df_long <- tidyr::pivot_longer(df_long, -Genus, names_to = "Sample", values_to = "Z")
  # keep sample order consistent with annotation
  df_long$Sample <- factor(df_long$Sample, levels = rownames(ann))
  
  # simple z-score palette roughly matching above
  zcols <- colorRampPalette(c("#073b4c", "#ffd166", "#ef476f"))(ncols)
  png(file.path(outdir, "heatmap_genus_top30_zscore.png"), width = 1600, height = 900, res = 140)
  ggplot(df_long, aes(x = Sample, y = Genus, fill = Z)) +
    geom_tile() +
    scale_fill_gradientn(colors = zcols) +
    labs(title = "Top 30 genera • row z-score of relative abundance", x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "right")
  dev.off()
}





##heatmap poption 3
# ensure output dir exists
outdir <- "exports"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 1) strip any leftover S4/sample_data class from the annotations
ann <- as.data.frame(ann)                    # drop S4 class
ann[] <- lapply(ann, function(x) {           # keep categorical as factors
  if (is.factor(x)) droplevels(x) else x
})

# 2) belt & suspenders: coerce the matrix and give explicit breaks
mat_z <- as.matrix(mat_z)
storage.mode(mat_z) <- "double"
mat_z[!is.finite(mat_z)] <- 0

ncols <- 101
cols_z <- colorRampPalette(c("#073b4c", "#ffd166", "#ef476f"))(ncols)
rng <- range(mat_z, na.rm = TRUE)
if (diff(rng) < .Machine$double.eps) rng <- rng + c(-1e-6, 1e-6)
brks <- seq(rng[1], rng[2], length.out = ncols + 1)

png(file.path(outdir, "heatmap_genus_top30_zscore.png"), width = 1600, height = 900, res = 140)
pheatmap::pheatmap(
  mat_z,
  color = cols_z,
  breaks = brks,
  annotation_col    = ann,
  annotation_colors = ann_colors,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  cluster_cols      = FALSE,
  clustering_method = "complete",
  clustering_distance_rows = "euclidean",
  main              = "Top 30 genera • row z-score of relative abundance"
)
dev.off()

utils::browseURL(normalizePath(file.path(outdir, "heatmap_genus_top30_zscore.png")))

# --- Make annotation plain & synced ---
outdir <- "exports"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ann_df <- as.data.frame(ann)                # drop S4 class
class(ann_df) <- "data.frame"               # hard reset class
ann_df$Zone <- droplevels(factor(as.character(ann_df$Zone)))
ann_df$Abbreviation <- droplevels(factor(as.character(ann_df$Abbreviation)))

stopifnot(identical(colnames(mat_z), rownames(ann_df)))  # columns == samples

# --- Rebuild palettes from *current* levels ---
zones <- levels(ann_df$Zone)
abbrs <- levels(ann_df$Abbreviation)

if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", repos="https://cloud.r-project.org")
make_big_palette <- function(n){ hues <- seq(15, 375, length.out = n + 1); grDevices::hcl(h=hues, l=65, c=100)[1:n] }

zone_cols <- setNames(RColorBrewer::brewer.pal(max(3, min(8, length(zones))), "Set2")[seq_along(zones)], zones)
abbr_cols <- setNames(make_big_palette(length(abbrs)), abbrs)

# --- Load ComplexHeatmap / circlize (install once if needed) ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
if (!requireNamespace("circlize", quietly = TRUE)) BiocManager::install("circlize", update = FALSE, ask = FALSE)

# --- Tiles color function (continuous for z-scores) ---
mat_z <- as.matrix(mat_z); storage.mode(mat_z) <- "double"; mat_z[!is.finite(mat_z)] <- 0
col_fun <- circlize::colorRamp2(
  c(min(mat_z, na.rm=TRUE), 0, max(mat_z, na.rm=TRUE)),
  c("#073b4c", "#ffd166", "#ef476f")
)

# --- Build annotation by columns (avoids the df= path that errored) ---
ha <- ComplexHeatmap::HeatmapAnnotation(
  Zone         = ann_df$Zone,
  Abbreviation = ann_df$Abbreviation,
  col          = list(Zone = zone_cols, Abbreviation = abbr_cols),
  which        = "column"
)

# --- Draw & save ---
png(file.path(outdir, "heatmap_genus_top30_zscore.png"), width = 1600, height = 900, res = 140)
print(ComplexHeatmap::Heatmap(
  mat_z,
  name              = "z",
  top_annotation    = ha,
  cluster_columns   = FALSE,     # keep your Zone→Abbreviation order
  show_column_names = FALSE,
  show_row_names    = TRUE,
  row_dend_side     = "left",
  column_title      = "Top 30 genera • row z-score of relative abundance",
  col               = col_fun
))
dev.off()

utils::browseURL(normalizePath(file.path(outdir, "heatmap_genus_top30_zscore.png")))







####Heatmap option 4
# install once if needed

# --- Install (once) ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (pkg in c("ComplexHeatmap","circlize")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

# --- Make sure inputs are plain + numeric ---
outdir <- "exports"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ann   <- as.data.frame(ann)                   # drop S4 class from sample_data
ann[] <- lapply(ann, function(x) if (is.factor(x)) droplevels(x) else x)

mat_z <- as.matrix(mat_z); storage.mode(mat_z) <- "double"
mat_z[!is.finite(mat_z)] <- 0                 # clean NAs/Infs

# --- Colors: continuous for tiles; discrete for annotations ---
col_fun <- circlize::colorRamp2(
  c(min(mat_z, na.rm = TRUE), 0, max(mat_z, na.rm = TRUE)),
  c("#073b4c", "#ffd166", "#ef476f")
)

ha <- ComplexHeatmap::HeatmapAnnotation(
  df  = ann[, c("Zone","Abbreviation"), drop = FALSE],
  col = list(Zone = zone_cols, Abbreviation = abbr_cols),
  which = "column"
)

# --- Draw & save ---
png(file.path(outdir, "heatmap_genus_top30_zscore.png"), width = 1600, height = 900, res = 140)
print(ComplexHeatmap::Heatmap(
  mat_z,
  name              = "z",
  top_annotation    = ha,
  cluster_columns   = FALSE,          # keep your Zone→Abbreviation order
  show_column_names = FALSE,
  show_row_names    = TRUE,
  row_dend_side     = "left",
  column_title      = "Top 30 genera • row z-score of relative abundance",
  col               = col_fun
))
dev.off()

utils::browseURL(normalizePath(file.path(outdir, "heatmap_genus_top30_zscore.png")))



# --- deps ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
for (pkg in c("ComplexHeatmap","circlize","RColorBrewer")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
}
outdir <- "exports"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- column annotation (no gradients; Zone & Abbreviation are discrete) ---
sdat <- as.data.frame(sample_data(ps_use))
stopifnot(all(c("Zone","Abbreviation") %in% names(sdat)))
ann_df <- sdat[, c("Zone","Abbreviation"), drop = FALSE]
ann_df$Zone <- droplevels(factor(as.character(ann_df$Zone)))
ann_df$Abbreviation <- droplevels(factor(as.character(ann_df$Abbreviation)))

zones <- levels(ann_df$Zone)
abbrs <- levels(ann_df$Abbreviation)
make_big_palette <- function(n){ hues <- seq(15, 375, length.out = n + 1); grDevices::hcl(h=hues, l=65, c=100)[1:n] }
zone_cols <- setNames(RColorBrewer::brewer.pal(max(3, min(8, length(zones))), "Set2")[seq_along(zones)], zones)
abbr_cols <- setNames(make_big_palette(length(abbrs)), abbrs)

# --- genus matrix (relative abundance), top 30 ---
ps_genus   <- tax_glom(ps_use, taxrank = "Genus", NArm = TRUE)
ps_genus_r <- transform_sample_counts(ps_genus, function(x) x / sum(x))
mat <- as.matrix(otu_table(ps_genus_r)); if (!taxa_are_rows(ps_genus_r)) mat <- t(mat)
g_names <- as.vector(tax_table(ps_genus_r)[, "Genus"]); g_names[is.na(g_names) | g_names==""] <- "Unclassified"
rownames(mat) <- g_names
mat <- rowsum(mat, group = rownames(mat))
top_rows <- head(order(rowMeans(mat), decreasing = TRUE), n = min(30, nrow(mat)))
mat_top <- mat[top_rows, rownames(ann_df), drop = FALSE]   # align to sample order

# --- log scale & draw ---
mat_log <- log10(mat_top * 1000 + 1)  # keeps absolute differences; spreads small values
storage.mode(mat_log) <- "double"; mat_log[!is.finite(mat_log)] <- 0

col_fun <- circlize::colorRamp2(
  c(min(mat_log, na.rm=TRUE), median(mat_log, na.rm=TRUE), max(mat_log, na.rm=TRUE)),
  c("#073b4c", "#ffd166", "#ef476f")
)
ha <- ComplexHeatmap::HeatmapAnnotation(
  Zone         = ann_df$Zone,
  Abbreviation = ann_df$Abbreviation,
  col          = list(Zone = zone_cols, Abbreviation = abbr_cols),
  which        = "column"
)
png(file.path(outdir, "heatmap_genus_top30_log.png"), width = 1600, height = 900, res = 140)
print(ComplexHeatmap::Heatmap(
  mat_log, name="log10(rel×1000+1)", top_annotation=ha,
  cluster_columns=FALSE, show_column_names=FALSE, show_row_names=TRUE,
  row_dend_side="left", column_title="Top 30 genera • log10(relative abundance × 1000 + 1)",
  col=col_fun
))
dev.off()
utils::browseURL(normalizePath(file.path(outdir, "heatmap_genus_top30_log.png")))



###Ordination (PCoA Bray–Curtis; color = Abbreviation, facet = Zone)

# deps
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos="https://cloud.r-project.org")
if (!requireNamespace("dplyr", quietly = TRUE))  install.packages("dplyr",  repos="https://cloud.r-project.org")
if (!requireNamespace("tidyr", quietly = TRUE))  install.packages("tidyr",  repos="https://cloud.r-project.org")
library(ggplot2); library(dplyr); library(tidyr)

# long table of relative abundance by Genus
ps_genus_r <- transform_sample_counts(ps_genus, function(x) x / sum(x))
df <- phyloseq::psmelt(ps_genus_r)                    # Sample, Abundance, Genus, Zone, Abbreviation, ...
df$Genus <- as.character(df$Genus); df$Genus[is.na(df$Genus) | df$Genus==""] <- "Unclassified"

# keep top 20 genera; collapse others
topG <- names(sort(tapply(df$Abundance, df$Genus, mean, na.rm=TRUE), decreasing=TRUE))[1:min(20, length(unique(df$Genus)))]
df$Genus2 <- ifelse(df$Genus %in% topG, df$Genus, "Other")

# aggregate to Sample × Genus2 (still relative; sums to 1 per sample)
bar_df <- df %>%
  dplyr::count(Sample, Zone, Abbreviation, Genus2, wt = Abundance, name = "RelAbund") %>%
  dplyr::group_by(Sample) %>% dplyr::mutate(RelAbund = RelAbund / sum(RelAbund)) %>% dplyr::ungroup()

# wide palette for genera
genus_levels <- unique(bar_df$Genus2)
genus_cols <- setNames(make_big_palette(length(genus_levels)), genus_levels)
genus_cols["Other"] <- "#cccccc"

# order samples by Zone then Abbreviation
bar_df <- bar_df %>%
  dplyr::mutate(Sample = factor(Sample, levels = rownames(ann_df)),
                Zone = droplevels(Zone), Abbreviation = droplevels(Abbreviation))

p <- ggplot(bar_df, aes(x = Sample, y = RelAbund, fill = Genus2)) +
  geom_col(width = 0.95) +
  facet_grid(. ~ Zone, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = genus_cols) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Relative abundance", fill = "Genus",
       title = "Genus composition per sample • faceted by Zone") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.spacing.x = unit(6, "pt"),
        legend.position = "right")

ggsave(file.path(outdir, "barplot_genus_by_sample_faceted_zone.png"), p, width = 16, height = 9, dpi = 140)
utils::browseURL(normalizePath(file.path(outdir, "barplot_genus_by_sample_faceted_zone.png")))

###nmds jaccard

if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan", repos="https://cloud.r-project.org")  # for distances/ellipses
library(ggplot2)

# Bray–Curtis on relative abundance
ps_rel <- transform_sample_counts(ps_use, function(x) x / sum(x))
ord_pcoa <- ordinate(ps_rel, method = "PCoA", distance = "bray")

# get % variance for axes (if available)
var_expl <- try({
  eig <- ord_pcoa$values$Relative_eig
  paste0(round(100*eig[1:2], 1), "%")
}, silent = TRUE)
xlab <- if (inherits(var_expl, "try-error")) "PCoA1" else paste0("PCoA1 (", var_expl[1], ")")
ylab <- if (inherits(var_expl, "try-error")) "PCoA2" else paste0("PCoA2 (", var_expl[2], ")")

df_pcoa <- plot_ordination(ps_rel, ord_pcoa, type = "samples", justDF = TRUE)
df_pcoa$Zone <- droplevels(factor(as.character(df_pcoa$Zone)))
df_pcoa$Abbreviation <- droplevels(factor(as.character(df_pcoa$Abbreviation)))

# reuse Abbreviation palette (categorical, no gradients)
if (!exists("abbr_cols")) {
  abbrs <- levels(df_pcoa$Abbreviation)
  abbr_cols <- setNames(make_big_palette(length(abbrs)), abbrs)
}

p_pcoa <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, color = Abbreviation)) +
  geom_point(size = 3, alpha = 0.9) +
  facet_wrap(~ Zone, scales = "free") +
  scale_color_manual(values = abbr_cols) +
  labs(title = "PCoA (Bray–Curtis) • color = Abbreviation • facet = Zone",
       x = xlab, y = ylab, color = "Abbreviation") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")
ggsave(file.path(outdir, "ordination_pcoa_bray_by_zone.png"), p_pcoa, width = 12, height = 7, dpi = 140)
utils::browseURL(normalizePath(file.path(outdir, "ordination_pcoa_bray_by_zone.png")))


##nmds jaccard
set.seed(1)
ord_nmds <- ordinate(ps_rel, method = "NMDS", distance = "jaccard", try = 200, k = 2)
df_nmds <- plot_ordination(ps_rel, ord_nmds, type = "samples", justDF = TRUE)
df_nmds$Zone <- droplevels(factor(as.character(df_nmds$Zone)))
df_nmds$Abbreviation <- droplevels(factor(as.character(df_nmds$Abbreviation)))

p_nmds <- ggplot(df_nmds, aes(x = NMDS1, y = NMDS2, color = Abbreviation)) +
  geom_point(size = 3, alpha = 0.9) +
  facet_wrap(~ Zone, scales = "free") +
  scale_color_manual(values = abbr_cols) +
  labs(title = "NMDS (Jaccard) • color = Abbreviation • facet = Zone",
       x = "NMDS1", y = "NMDS2", color = "Abbreviation") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "ordination_nmds_jaccard_by_zone.png"), p_nmds, width = 12, height = 7, dpi = 140)
utils::browseURL(normalizePath(file.path(outdir, "ordination_nmds_jaccard_by_zone.png")))




# Core deps
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")
for (pkg in c("phyloseq","vegan","ggplot2","dplyr","tidyr","RColorBrewer","scales"))
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)

library(phyloseq); library(ggplot2); library(dplyr); library(tidyr)

outdir <- "exports"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Palettes: Abbreviation gets distinct categorical colors; Zone also categorical
make_big_palette <- function(n){
  hues <- seq(15, 375, length.out = n + 1)
  grDevices::hcl(h=hues, l=65, c=100)[1:n]
}
sdat <- as.data.frame(sample_data(ps_use))
stopifnot(all(c("Zone","Abbreviation") %in% names(sdat)))
sdat$Zone         <- droplevels(factor(as.character(sdat$Zone)))
sdat$Abbreviation <- droplevels(factor(as.character(sdat$Abbreviation)))
abbr_levels <- levels(sdat$Abbreviation)
zone_levels <- levels(sdat$Zone)
abbr_cols <- setNames(make_big_palette(length(abbr_levels)), abbr_levels)
zone_cols <- setNames(RColorBrewer::brewer.pal(max(3, min(8, length(zone_levels))), "Set2")[seq_along(zone_levels)], zone_levels)

# Keep a sample order (Zone -> Abbreviation) for consistent plotting
sample_order <- rownames(sdat[order(sdat$Zone, sdat$Abbreviation), , drop=FALSE])



# Core deps
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")
for (pkg in c("phyloseq","vegan","ggplot2","dplyr","tidyr","RColorBrewer","scales"))
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)

library(phyloseq); library(ggplot2); library(dplyr); library(tidyr)

outdir <- "exports"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Palettes: Abbreviation gets distinct categorical colors; Zone also categorical
make_big_palette <- function(n){
  hues <- seq(15, 375, length.out = n + 1)
  grDevices::hcl(h=hues, l=65, c=100)[1:n]
}
sdat <- as.data.frame(sample_data(ps_use))
stopifnot(all(c("Zone","Abbreviation") %in% names(sdat)))
sdat$Zone         <- droplevels(factor(as.character(sdat$Zone)))
sdat$Abbreviation <- droplevels(factor(as.character(sdat$Abbreviation)))
abbr_levels <- levels(sdat$Abbreviation)
zone_levels <- levels(sdat$Zone)
abbr_cols <- setNames(make_big_palette(length(abbr_levels)), abbr_levels)
zone_cols <- setNames(RColorBrewer::brewer.pal(max(3, min(8, length(zone_levels))), "Set2")[seq_along(zone_levels)], zone_levels)

# Keep a sample order (Zone -> Abbreviation) for consistent plotting
sample_order <- rownames(sdat[order(sdat$Zone, sdat$Abbreviation), , drop=FALSE])


# (Optional) rarefy BEFORE alpha if you want equal depth comparisons
DO_RAREFY <- FALSE; RARE_DEPTH <- NA
ps_alpha <- if (DO_RAREFY) {
  set.seed(1)
  rarefy_even_depth(ps_use, sample.size = if (is.na(RARE_DEPTH)) min(sample_sums(ps_use)) else RARE_DEPTH,
                    rngseed = 1, replace = FALSE, verbose = FALSE)
} else ps_use

alpha_df <- estimate_richness(ps_alpha, measures = c("Shannon","Simpson","Observed","Chao1")) %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(sdat %>% tibble::rownames_to_column("Sample"), by = "Sample") %>%
  mutate(Sample = factor(Sample, levels = sample_order))

plot_alpha <- function(metric){
  ggplot(alpha_df, aes(x = Zone, y = .data[[metric]], color = Abbreviation)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.85) +
    scale_color_manual(values = abbr_cols, drop = FALSE) +
    labs(title = paste0("Alpha diversity • ", metric),
         x = "Zone", y = metric, color = "Abbreviation") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
}

for (m in c("Shannon","Simpson","Observed","Chao1")) {
  p <- plot_alpha(m)
  ggsave(file.path(outdir, paste0("alpha_", tolower(m), ".png")), p, width=9, height=6, dpi=140)
}


# (Optional) rarefy BEFORE alpha if you want equal depth comparisons
DO_RAREFY <- FALSE; RARE_DEPTH <- NA
ps_alpha <- if (DO_RAREFY) {
  set.seed(1)
  rarefy_even_depth(ps_use, sample.size = if (is.na(RARE_DEPTH)) min(sample_sums(ps_use)) else RARE_DEPTH,
                    rngseed = 1, replace = FALSE, verbose = FALSE)
} else ps_use

alpha_df <- estimate_richness(ps_alpha, measures = c("Shannon","Simpson","Observed","Chao1")) %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(sdat %>% tibble::rownames_to_column("Sample"), by = "Sample") %>%
  mutate(Sample = factor(Sample, levels = sample_order))

plot_alpha <- function(metric){
  ggplot(alpha_df, aes(x = Zone, y = .data[[metric]], color = Abbreviation)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.85) +
    scale_color_manual(values = abbr_cols, drop = FALSE) +
    labs(title = paste0("Alpha diversity • ", metric),
         x = "Zone", y = metric, color = "Abbreviation") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
}

for (m in c("Shannon","Simpson","Observed","Chao1")) {
  p <- plot_alpha(m)
  ggsave(file.path(outdir, paste0("alpha_", tolower(m), ".png")), p, width=9, height=6, dpi=140)
}

# Use relative abundance for Bray–Curtis
ps_rel <- transform_sample_counts(ps_use, function(x) x / sum(x))
ord_pcoa <- ordinate(ps_rel, method = "PCoA", distance = "bray")

# % variance labels if available
eig <- ord_pcoa$values$Relative_eig
xlab <- if (!is.null(eig)) paste0("PCoA1 (", round(100*eig[1],1), "%)") else "PCoA1"
ylab <- if (!is.null(eig)) paste0("PCoA2 (", round(100*eig[2],1), "%)") else "PCoA2"

df_pcoa <- plot_ordination(ps_rel, ord_pcoa, type = "samples", justDF = TRUE) %>%
  mutate(Zone = droplevels(factor(as.character(Zone))),
         Abbreviation = droplevels(factor(as.character(Abbreviation))),
         Sample = factor(Sample, levels = sample_order))

p_pcoa <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, color = Abbreviation)) +
  geom_point(size = 3, alpha = 0.9) +
  facet_wrap(~ Zone, scales = "free") +
  scale_color_manual(values = abbr_cols, drop = FALSE) +
  labs(title = "PCoA • Bray–Curtis", x = xlab, y = ylab, color = "Abbreviation") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "beta_pcoa_bray_by_zone.png"), p_pcoa, width = 12, height = 7, dpi = 140)


####if (!requireNamespace("ANCOMBC", quietly = TRUE))
  BiocManager::install("ANCOMBC", update = FALSE, ask = FALSE)
if (!requireNamespace("ANCOMBC2", quietly = TRUE))
  BiocManager::install("ANCOMBC2", update = FALSE, ask = FALSE)
library(ANCOMBC2)

# ---- ASV-level test (no aggregation) ----
res_asv <- ANCOMBC2::ancombc2(
  data            = ps_use,
  tax_level       = NULL,              # NULL = use ASVs
  fix_formula     = "Zone",            # model ~ Zone
  rand_formula    = NULL,
  p_adj_method    = "holm",
  prv_cut         = 0.10,              # prevalence filter
  lib_cut         = 1000,              # library size filter
  s0_perc         = 0.05,
  group           = "Zone",            # for pairwise/contrast
  struc_zero      = TRUE,
  neg_lb          = TRUE,
  alpha           = 0.05,
  global          = TRUE
)
# Significant ASVs
da_asv <- res_asv$res %>% dplyr::filter(q_val < 0.05)
write.csv(res_asv$res, file.path(outdir, "DA_ANCOMBC2_ASV_full.csv"), row.names = FALSE)
write.csv(da_asv,       file.path(outdir, "DA_ANCOMBC2_ASV_sig_q<0.05.csv"), row.names = FALSE)

# ---- Genus-level test (aggregates to Genus) ----
# (Requires taxonomy with a "Genus" column)
if ("Genus" %in% colnames(tax_table(ps_use))) {
  res_genus <- ANCOMBC2::ancombc2(
    data            = ps_use,
    tax_level       = "Genus",
    fix_formula     = "Zone",
    rand_formula    = NULL,
    p_adj_method    = "holm",
    prv_cut         = 0.10,
    lib_cut         = 1000,
    s0_perc         = 0.05,
    group           = "Zone",
    struc_zero      = TRUE,
    neg_lb          = TRUE,
    alpha           = 0.05,
    global          = TRUE
  )
  da_genus <- res_genus$res %>% dplyr::filter(q_val < 0.05)
  write.csv(res_genus$res, file.path(outdir, "DA_ANCOMBC2_Genus_full.csv"), row.names = FALSE)
  write.csv(da_genus,       file.path(outdir, "DA_ANCOMBC2_Genus_sig_q<0.05.csv"), row.names = FALSE)
}

# Quick volcano-style plot for genus (if present)
if (exists("res_genus")) {
  vol <- res_genus$res %>% mutate(sig = q_val < 0.05)
  p_vol <- ggplot(vol, aes(x = lfc, y = -log10(p_val), color = sig)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("FALSE"="#999999","TRUE"="#d62728")) +
    labs(title = "ANCOM-BC2 (Genus) • Zone contrasts", x = "log fold-change", y = "-log10(p)") +
    theme_minimal(base_size = 12) + theme(legend.position="none")
  ggsave(file.path(outdir, "DA_ANCOMBC2_genus_volcano.png"), p_vol, width=8, height=6, dpi=140)
}

###ANCOM-BC2
if (!requireNamespace("ANCOMBC", quietly = TRUE))
  BiocManager::install("ANCOMBC", update = FALSE, ask = FALSE)
if (!requireNamespace("ANCOMBC2", quietly = TRUE))
  BiocManager::install("ANCOMBC2", update = FALSE, ask = FALSE)
library(ANCOMBC2)

# ---- ASV-level test (no aggregation) ----
res_asv <- ANCOMBC2::ancombc2(
  data            = ps_use,
  tax_level       = NULL,              # NULL = use ASVs
  fix_formula     = "Zone",            # model ~ Zone
  rand_formula    = NULL,
  p_adj_method    = "holm",
  prv_cut         = 0.10,              # prevalence filter
  lib_cut         = 1000,              # library size filter
  s0_perc         = 0.05,
  group           = "Zone",            # for pairwise/contrast
  struc_zero      = TRUE,
  neg_lb          = TRUE,
  alpha           = 0.05,
  global          = TRUE
)
# Significant ASVs
da_asv <- res_asv$res %>% dplyr::filter(q_val < 0.05)
write.csv(res_asv$res, file.path(outdir, "DA_ANCOMBC2_ASV_full.csv"), row.names = FALSE)
write.csv(da_asv,       file.path(outdir, "DA_ANCOMBC2_ASV_sig_q<0.05.csv"), row.names = FALSE)

# ---- Genus-level test (aggregates to Genus) ----
# (Requires taxonomy with a "Genus" column)
if ("Genus" %in% colnames(tax_table(ps_use))) {
  res_genus <- ANCOMBC2::ancombc2(
    data            = ps_use,
    tax_level       = "Genus",
    fix_formula     = "Zone",
    rand_formula    = NULL,
    p_adj_method    = "holm",
    prv_cut         = 0.10,
    lib_cut         = 1000,
    s0_perc         = 0.05,
    group           = "Zone",
    struc_zero      = TRUE,
    neg_lb          = TRUE,
    alpha           = 0.05,
    global          = TRUE
  )
  da_genus <- res_genus$res %>% dplyr::filter(q_val < 0.05)
  write.csv(res_genus$res, file.path(outdir, "DA_ANCOMBC2_Genus_full.csv"), row.names = FALSE)
  write.csv(da_genus,       file.path(outdir, "DA_ANCOMBC2_Genus_sig_q<0.05.csv"), row.names = FALSE)
}

# Quick volcano-style plot for genus (if present)
if (exists("res_genus")) {
  vol <- res_genus$res %>% mutate(sig = q_val < 0.05)
  p_vol <- ggplot(vol, aes(x = lfc, y = -log10(p_val), color = sig)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("FALSE"="#999999","TRUE"="#d62728")) +
    labs(title = "ANCOM-BC2 (Genus) • Zone contrasts", x = "log fold-change", y = "-log10(p)") +
    theme_minimal(base_size = 12) + theme(legend.position="none")
  ggsave(file.path(outdir, "DA_ANCOMBC2_genus_volcano.png"), p_vol, width=8, height=6, dpi=140)
}



###PHYLOGENETIC TREE
###



for (pkg in c("DECIPHER","phangorn"))
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
library(DECIPHER); library(phangorn)

# --- gather sequences ---
txdf <- as.data.frame(tax_table(ps_use))
seqs <- NULL
# If taxa names look like DNA strings, use them
tn <- taxa_names(ps_use)
if (all(grepl("^[ACGTN]+$", tn))) seqs <- tn
# Else try taxonomy column "Sequence"
if (is.null(seqs) && "Sequence" %in% names(txdf)) seqs <- txdf$Sequence
stopifnot(!is.null(seqs))
names(seqs) <- taxa_names(ps_use)

# --- multiple alignment & tree ---
dna <- Biostrings::DNAStringSet(seqs)
names(dna) <- taxa_names(ps_use)
alignment <- DECIPHER::AlignSeqs(dna, processors = max(1, parallel::detectCores()-1))

phydat <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
dm <- phangorn::dist.ml(phydat)
treeNJ <- phangorn::NJ(dm)
fit <- phangorn::pml(treeNJ, phydat)
fitOpt <- phangorn::optim.pml(fit, model="GTR", optInv=TRUE, optGamma=TRUE, control = phangorn::pml.control(trace = 0))
tree <- fitOpt$tree

# Add tree to phyloseq
ps_tree <- phyloseq(otu_table(ps_use), sample_data(ps_use), tax_table(ps_use), phy_tree(tree))

# UniFrac ordination (weighted + unweighted)
ord_w  <- ordinate(ps_tree, method = "PCoA", distance = "wunifrac")
ord_uw <- ordinate(ps_tree, method = "PCoA", distance = "unifrac")

plot_uni <- function(ord_obj, title_txt){
  df <- plot_ordination(ps_tree, ord_obj, type="samples", justDF=TRUE) %>%
    mutate(Zone = droplevels(factor(as.character(Zone))),
           Abbreviation = droplevels(factor(as.character(Abbreviation))))
  eig <- ord_obj$values$Relative_eig
  xlab <- if (!is.null(eig)) paste0("PCoA1 (", round(100*eig[1],1), "%)") else "PCoA1"
  ylab <- if (!is.null(eig)) paste0("PCoA2 (", round(100*eig[2],1), "%)") else "PCoA2"
  ggplot(df, aes(Axis.1, Axis.2, color = Abbreviation)) +
    geom_point(size=3, alpha=0.9) +
    facet_wrap(~ Zone, scales = "free") +
    scale_color_manual(values = abbr_cols, drop = FALSE) +
    labs(title = title_txt, x = xlab, y = ylab, color = "Abbreviation") +
    theme_minimal(base_size = 12)
}

p_w  <- plot_uni(ord_w,  "PCoA • Weighted UniFrac")
p_uw <- plot_uni(ord_uw, "PCoA • Unweighted UniFrac")
ggsave(file.path(outdir, "beta_pcoa_wunifrac_by_zone.png"),  p_w,  width=12, height=7, dpi=140)
ggsave(file.path(outdir, "beta_pcoa_unifrac_by_zone.png"),   p_uw, width=12, height=7, dpi=140)


####
####
####


# Install what’s missing
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

pkgs_cran <- c("ape", "phangorn")
for(p in pkgs_cran)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

for(p in c("DECIPHER","Biostrings"))
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, update = FALSE, ask = FALSE)

library(DECIPHER); library(Biostrings); library(ape); library(phangorn)

## --- Collect ASV sequences (same logic as before) ---
txdf <- as.data.frame(tax_table(ps_use))
seqs <- NULL
tn <- taxa_names(ps_use)
if (all(grepl("^[ACGTN]+$", tn))) seqs <- tn
if (is.null(seqs) && "Sequence" %in% names(txdf)) seqs <- txdf$Sequence
stopifnot(!is.null(seqs))
names(seqs) <- taxa_names(ps_use)

## --- Align & infer tree with phangorn ---
dna <- Biostrings::DNAStringSet(seqs)
names(dna) <- taxa_names(ps_use)
alignment <- DECIPHER::AlignSeqs(dna, processors = max(1, parallel::detectCores()-1))

phydat <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
dm     <- phangorn::dist.ml(phydat)
treeNJ <- phangorn::NJ(dm)
fit    <- phangorn::pml(treeNJ, phydat)
fitOpt <- phangorn::optim.pml(fit, model="GTR", optInv=TRUE, optGamma=TRUE,
                              control = phangorn::pml.control(trace = 0))
tree   <- fitOpt$tree

ps_tree <- phyloseq(otu_table(ps_use), sample_data(ps_use), tax_table(ps_use), phy_tree(tree))

# 0) Optional: silence the 'phylo' cache message
if ("package:tidytree" %in% search()) detach("package:tidytree", unload = TRUE)

# 1) Ensure the tree is rooted (midpoint) and has valid edge lengths
if (!requireNamespace("ape", quietly = TRUE))
  install.packages("ape", repos = "https://cloud.r-project.org")
if (!requireNamespace("phangorn", quietly = TRUE))
  install.packages("phangorn", repos = "https://cloud.r-project.org")

library(ape); library(phangorn)

# Root at the midpoint if not already rooted
if (!ape::is.rooted(tree)) {
  tree <- phangorn::midpoint(tree)
}

# Make sure all branch lengths are positive (UniFrac needs this)
if (!is.null(tree$edge.length)) {
  tree$edge.length[tree$edge.length <= 0 | !is.finite(tree$edge.length)] <- 1e-8
}

# 2) Rebuild the phyloseq with the rooted tree
ps_tree <- phyloseq(otu_table(ps_use),
                    sample_data(ps_use),
                    tax_table(ps_use),
                    phy_tree(tree))

# 3) Recompute UniFrac ordinations (no more random root warning)
ord_w  <- ordinate(ps_tree, method = "PCoA", distance = "wunifrac")
ord_uw <- ordinate(ps_tree, method = "PCoA", distance = "unifrac")

# 4) Plot (color = Abbreviation, faceted by Zone — no gradients)
plot_uni <- function(ord_obj, title_txt){
  df <- plot_ordination(ps_tree, ord_obj, type = "samples", justDF = TRUE)
  ggplot(df, aes(Axis.1, Axis.2, color = Abbreviation)) +
    geom_point(size = 3) +
    facet_wrap(~ Zone, scales = "free") +
    scale_color_manual(
      values = setNames(grDevices::hcl(seq(15, 375, length.out = length(levels(sample_data(ps_tree)$Abbreviation)) + 1)[-1],
                                       l = 65, c = 100),
                        levels(sample_data(ps_tree)$Abbreviation))
    ) +
    labs(title = title_txt, x = "PCoA1", y = "PCoA2", color = "Abbreviation") +
    theme_minimal(base_size = 12)
}

dir.create("exports", showWarnings = FALSE)
ggsave("exports/beta_pcoa_wunifrac_by_zone.png",  plot_uni(ord_w,  "PCoA • Weighted UniFrac"),  width = 12, height = 7, dpi = 140)
ggsave("exports/beta_pcoa_unifrac_by_zone.png",   plot_uni(ord_uw, "PCoA • Unweighted UniFrac"), width = 12, height = 7, dpi = 140)





## UniFrac ordinations (saves figures in exports/)
ord_w  <- ordinate(ps_tree, method="PCoA", distance="wunifrac")
ord_uw <- ordinate(ps_tree, method="PCoA", distance="unifrac")

plot_uni <- function(ord_obj, title_txt){
  df <- plot_ordination(ps_tree, ord_obj, type="samples", justDF=TRUE)
  ggplot(df, aes(Axis.1, Axis.2, color = Abbreviation)) +
    geom_point(size=3) +
    facet_wrap(~ Zone, scales = "free") +
    scale_color_manual(values = setNames(grDevices::hcl(seq(15,375,length.out=length(levels(sample_data(ps_tree)$Abbreviation))+1)[-1],
                                                        l=65,c=100),
                                         levels(sample_data(ps_tree)$Abbreviation))) +
    labs(title = title_txt, x = "PCoA1", y = "PCoA2", color = "Abbreviation") +
    theme_minimal(base_size = 12)
}
dir.create("exports", showWarnings = FALSE)
ggsave("exports/beta_pcoa_wunifrac_by_zone.png",  plot_uni(ord_w,  "PCoA • Weighted UniFrac"),  width=12, height=7, dpi=140)
ggsave("exports/beta_pcoa_unifrac_by_zone.png",   plot_uni(ord_uw, "PCoA • Unweighted UniFrac"), width=12, height=7, dpi=140)




#####
#####
#####


















# ---------- Alpha diversity (x = Zone; color = Abbreviation) ----------
alpha_df <- estimate_richness(ps_use, measures = c("Shannon","Simpson","Observed")) %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(meta_ps %>% tibble::rownames_to_column("Sample"), by = "Sample")
write.csv(alpha_df, file.path(outdir, "alpha_diversity_byZoneAbbrev.csv"), row.names = FALSE)

for(metric in c("Shannon","Simpson","Observed")){
  gp <- ggplot(alpha_df, aes(x = .data[[ZONE_COL]], y = .data[[metric]], color = .data[[ABBREV_COL]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.95) +
    labs(title = paste("Alpha diversity:", metric), x = "Zone", y = metric, color = "Individual") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  ggsave(file.path(outdir, paste0("alpha_", tolower(metric), "_byZone_colorAbbrev.png")), gp, width = 8, height = 5, dpi = 300)
}

# ---------- Beta diversity & ordinations (color = Abbreviation; shape = Zone) ----------
dist_bray <- phyloseq::distance(ps_use, method="bray")
set.seed(1)
ord_nmds <- ordinate(ps_use, method="NMDS", distance=dist_bray, trymax=100)
ord_pcoa <- ordinate(ps_use, method="PCoA", distance=dist_bray)

p_nmds <- plot_ordination(ps_use, ord_nmds, type="samples", color=ABBREV_COL, shape=ZONE_COL) +
  geom_point(size=3, alpha=0.95) +
  labs(title="NMDS (Bray-Curtis): color=Individual, shape=Zone", color="Individual", shape="Zone")
ggsave(file.path(outdir, "ordination_NMDS_bray_colorAbbrev_shapeZone.png"), p_nmds, width=7, height=5, dpi=300)

p_pcoa <- plot_ordination(ps_use, ord_pcoa, type="samples", color=ABBREV_COL, shape=ZONE_COL) +
  geom_point(size=3, alpha=0.95) +
  labs(title="PCoA (Bray-Curtis): color=Individual, shape=Zone", color="Individual", shape="Zone")
ggsave(file.path(outdir, "ordination_PCoA_bray_colorAbbrev_shapeZone.png"), p_pcoa, width=7, height=5, dpi=300)

# Save ordination coordinates
ord_scores <- as.data.frame(scores(ord_nmds, display="sites")) %>% tibble::rownames_to_column("Sample")
ord_scores <- left_join(ord_scores, meta_ps %>% tibble::rownames_to_column("Sample"), by="Sample")
write.csv(ord_scores, file.path(outdir, "ordination_NMDS_scores_withZoneAbbrev.csv"), row.names = FALSE)

pcoa_scores <- as.data.frame(ord_pcoa$vectors) %>% tibble::rownames_to_column("Sample")
pcoa_scores <- left_join(pcoa_scores, meta_ps %>% tibble::rownames_to_column("Sample"), by="Sample")
write.csv(pcoa_scores, file.path(outdir, "ordination_PCoA_scores_withZoneAbbrev.csv"), row.names = FALSE)

# PERMANOVA by Zone (primary); also a version stratified by Abbreviation (if repeated measures)
set.seed(1)
adonis_zone <- vegan::adonis2(as.matrix(dist_bray) ~ meta_ps[[ZONE_COL]], permutations=999)
capture.output(adonis_zone, file=file.path(outdir, "PERMANOVA_Zone.txt"))

if(length(unique(meta_ps[[ABBREV_COL]])) < nrow(meta_ps)){
  set.seed(1)
  adonis_zone_block <- vegan::adonis2(as.matrix(dist_bray) ~ meta_ps[[ZONE_COL]],
                                      permutations = how(nperm=999, strata=meta_ps[[ABBREV_COL]]))
  capture.output(adonis_zone_block, file=file.path(outdir, "PERMANOVA_Zone_stratifiedByAbbreviation.txt"))
}

# ---------- Taxonomic barplots (if taxonomy present) ----------
if(!all(is.na(tax_table(ps_use)[,1]))){
  ps_genus <- tax_glom(ps_use, taxrank="Genus", NArm=TRUE)
  ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x/sum(x))
  topN <- 12
  top_genera <- names(sort(taxa_sums(ps_genus_rel), decreasing=TRUE))[1:min(topN, ntaxa(ps_genus_rel))]
  df_genus <- psmelt(prune_taxa(top_genera, ps_genus_rel))
  df_genus[[ZONE_COL]] <- sample_data(ps_genus_rel)[[ZONE_COL]][ match(df_genus$Sample, rownames(sample_data(ps_genus_rel))) ]
  
  gp_genus <- ggplot(df_genus, aes(x=Sample, y=Abundance, fill=Genus)) +
    geom_col() + facet_wrap(as.formula(paste("~", ZONE_COL)), scales="free_x") +
    labs(title="Relative Abundance by Genus (faceted by Zone)", y="Relative Abundance") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
  ggsave(file.path(outdir, "barplot_genus_relabund_facetedByZone.png"), gp_genus, width=12, height=6, dpi=300)
} else {
  ps_rel <- transform_sample_counts(ps_use, function(x) x/sum(x))
  rel_long <- psmelt(ps_rel) %>%
    dplyr::group_by(OTU) %>% summarize(mean_rel = mean(Abundance), .groups="drop") %>%
    arrange(desc(mean_rel)) %>% slice(1:20)
  write.csv(rel_long, file.path(outdir, "top20_ASVs_relative.csv"), row.names = FALSE)
}

message("Done. See the 'exports/' folder for results.")


####ANALYSIS BASED ON LOCATION 
####
## =========================
## Location-averaged Diversity (Alpha & Beta)
## =========================
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
})

# ---- CONFIG: set your column names here if needed ----
LOCATION_COL   <- if (exists("LOCATION_COL")) LOCATION_COL else "Location"
ZONE_COL       <- if (exists("ZONE_COL"))     ZONE_COL     else "Zone"
ABBREV_COL     <- if (exists("ABBREV_COL"))   ABBREV_COL   else "Abbreviation"

stopifnot(inherits(ps_use, "phyloseq"))
sd_ps <- as(sample_data(ps_use), "data.frame")
stopifnot(all(c(LOCATION_COL, ZONE_COL, ABBREV_COL) %in% names(sd_ps)))

# (Optional) enforce that we have all 29 samples
if (nsamples(ps_use) != 29) {
  warning(sprintf("Found %d samples in ps_use (expected 29). Proceeding with all available.",
                  nsamples(ps_use)))
}

# Make sure groups are factors with stable levels
sd_ps[[LOCATION_COL]] <- factor(sd_ps[[LOCATION_COL]])
sd_ps[[ZONE_COL]]     <- factor(sd_ps[[ZONE_COL]])
sd_ps[[ABBREV_COL]]   <- factor(sd_ps[[ABBREV_COL]])
sample_data(ps_use)   <- sample_data(sd_ps)

# Reorder Locations by Zone first, then alphabetically within Zone
loc_by_zone <- sd_ps %>%
  count(!!sym(LOCATION_COL), !!sym(ZONE_COL), sort = FALSE) %>%
  group_by(!!sym(LOCATION_COL)) %>%
  slice_max(n, with_ties = FALSE) %>%        # dominant Zone for that Location
  ungroup() %>%
  arrange(!!sym(ZONE_COL), !!sym(LOCATION_COL))

sd_ps[[LOCATION_COL]] <- factor(sd_ps[[LOCATION_COL]],
                                levels = unique(loc_by_zone[[LOCATION_COL]]))
sample_data(ps_use) <- sample_data(sd_ps)

# Palettes: distinct colors for Abbreviation (individuals), shapes for Zone (no gradient)
make_big_palette <- function(n){
  hues <- seq(15, 375, length.out = n + 1)
  hcl(h = hues[-length(hues)], l = 65, c = 100)
}
abbr_levels <- levels(sd_ps[[ABBREV_COL]])
zone_levels <- levels(sd_ps[[ZONE_COL]])
col_abbr    <- setNames(make_big_palette(length(abbr_levels)), abbr_levels)

# Up to 6 distinct shapes. Add more if you have more Zones.
shape_vals  <- c(16, 17, 15, 3, 8, 7)[seq_along(zone_levels)]
names(shape_vals) <- zone_levels

# Output dir
outdir <- "exports"
dir.create(outdir, showWarnings = FALSE)

## =========================
## ALPHA DIVERSITY
## =========================
# Compute per-sample alpha metrics
alpha_df <- estimate_richness(
  ps_use, measures = c("Observed","Shannon","Simpson","Chao1")
) %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(sd_ps %>% tibble::rownames_to_column("SampleID"),
            by = "SampleID") %>%
  relocate(all_of(c("SampleID", LOCATION_COL, ZONE_COL, ABBREV_COL)), .before = 1)

# Per-Location averages (mean ± sd)
alpha_loc_avg <- alpha_df %>%
  group_by(!!sym(LOCATION_COL)) %>%
  summarize(
    n = dplyr::n(),
    Observed_mean = mean(Observed, na.rm = TRUE),
    Observed_sd   = sd(Observed, na.rm = TRUE),
    Shannon_mean  = mean(Shannon, na.rm = TRUE),
    Shannon_sd    = sd(Shannon, na.rm = TRUE),
    Simpson_mean  = mean(Simpson, na.rm = TRUE),
    Simpson_sd    = sd(Simpson, na.rm = TRUE),
    Chao1_mean    = mean(Chao1, na.rm = TRUE),
    Chao1_sd      = sd(Chao1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # keep the Location ordering we defined
  mutate(!!LOCATION_COL := factor(.data[[LOCATION_COL]],
                                  levels = levels(sd_ps[[LOCATION_COL]])))

# Save tables
write.csv(alpha_df,      file.path(outdir, "alpha_per_sample.csv"), row.names = FALSE)
write.csv(alpha_loc_avg, file.path(outdir, "alpha_location_averages.csv"), row.names = FALSE)

# Plot helper for alpha diversity
plot_alpha_metric <- function(metric, ylab_txt){
  ggplot(alpha_df,
         aes(x = .data[[LOCATION_COL]],
             y = .data[[metric]],
             color = .data[[ABBREV_COL]],
             shape = .data[[ZONE_COL]])) +
    # individual points (identify each sample)
    geom_point(position = position_jitter(width = 0.12, height = 0),
               size = 2.6, alpha = 0.95) +
    # optional labels: use sample Abbreviation; switch to SampleID if you prefer
    ggrepel::geom_text_repel(aes(label = .data[[ABBREV_COL]]),
                             size = 2.4, max.overlaps = 20,
                             position = position_jitter(width = 0.12, height = 0),
                             segment.color = "grey80", alpha = 0.8, show.legend = FALSE) +
    # Location mean ± SE line/point
    stat_summary(aes(group = 1),
                 fun = mean, geom = "line", color = "black", size = 0.5) +
    stat_summary(fun.data = mean_se, geom = "pointrange",
                 color = "black", size = 0.4) +
    scale_color_manual(values = col_abbr, name = "Individual") +
    scale_shape_manual(values = shape_vals, name = "Zone") +
    labs(x = "Location (ordered by Zone)", y = ylab_txt,
         title = paste("Alpha diversity:", metric),
         subtitle = "Points = samples (colored by individual, shaped by Zone); black = Location mean ± SE") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

# Save alpha plots for key metrics
ggsave(file.path(outdir, "alpha_Shannon_by_location.png"),
       plot_alpha_metric("Shannon", "Shannon index"),
       width = 12, height = 6, dpi = 140)

ggsave(file.path(outdir, "alpha_Chao1_by_location.png"),
       plot_alpha_metric("Chao1", "Chao1 richness"),
       width = 12, height = 6, dpi = 140)

ggsave(file.path(outdir, "alpha_Observed_by_location.png"),
       plot_alpha_metric("Observed", "Observed ASVs"),
       width = 12, height = 6, dpi = 140)

## =========================
## BETA DIVERSITY (Bray–Curtis) with Location centroids
## =========================
# Use relative abundance for composition distances

## =========================
## Diversity & DA master script
## =========================

## ---- Package setup ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

if (!requireNamespace("phyloseq", quietly = TRUE))
  BiocManager::install("phyloseq", update = FALSE, ask = FALSE)

for (p in c("ggplot2","ggrepel","dplyr","tidyr","vegan","readr")) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ANCOMBC", quietly = TRUE))
  BiocManager::install("ANCOMBC", update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(phyloseq); library(dplyr); library(tidyr)
  library(ggplot2); library(ggrepel); library(vegan); library(ANCOMBC)
})

## ---- Config (edit here if your column names differ) ----
LOCATION_COL <- "Location"
ZONE_COL     <- "Zone"
ABBREV_COL   <- "Abbreviation"

stopifnot(inherits(ps_use, "phyloseq"))
sd_ps <- as(sample_data(ps_use), "data.frame")
need <- c(LOCATION_COL, ZONE_COL, ABBREV_COL)
stopifnot(all(need %in% names(sd_ps)))

# Optional sanity check for 29 samples
if (nsamples(ps_use) != 29) {
  warning(sprintf("Found %d samples in ps_use (expected 29). Proceeding with all available.",
                  nsamples(ps_use)))
}

# Factors with stable levels
sd_ps[[LOCATION_COL]] <- factor(sd_ps[[LOCATION_COL]])
sd_ps[[ZONE_COL]]     <- factor(sd_ps[[ZONE_COL]])
sd_ps[[ABBREV_COL]]   <- factor(sd_ps[[ABBREV_COL]])
sample_data(ps_use)   <- sample_data(sd_ps)

# Order Locations by (Zone → Location) using a dominant Zone per Location
loc_by_zone <- sd_ps %>%
  count(!!rlang::sym(LOCATION_COL), !!rlang::sym(ZONE_COL), name = "n", sort = FALSE) %>%
  group_by(!!rlang::sym(LOCATION_COL)) %>%
  slice_max(n, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(!!rlang::sym(ZONE_COL), !!rlang::sym(LOCATION_COL))

sd_ps[[LOCATION_COL]] <- factor(sd_ps[[LOCATION_COL]], levels = unique(loc_by_zone[[LOCATION_COL]]))
sample_data(ps_use)   <- sample_data(sd_ps)

# Palettes: distinct colors for individuals, shapes for Zone
make_big_palette <- function(n){
  hues <- seq(15, 375, length.out = n + 1)
  grDevices::hcl(h = hues[-length(hues)], l = 65, c = 100)
}
abbr_levels <- levels(sd_ps[[ABBREV_COL]])
zone_levels <- levels(sd_ps[[ZONE_COL]])
col_abbr    <- setNames(make_big_palette(length(abbr_levels)), abbr_levels)
shape_vals  <- setNames(rep(c(16,17,15,3,8,7,4,18,0,1), length.out = length(zone_levels)), zone_levels)

# Output dir
outdir <- "exports"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## =========================
## ALPHA DIVERSITY
## =========================
alpha_df <- estimate_richness(ps_use, measures = c("Observed","Shannon","Simpson","Chao1")) %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(sd_ps %>% tibble::rownames_to_column("SampleID"), by = "SampleID") %>%
  relocate(all_of(c("SampleID", LOCATION_COL, ZONE_COL, ABBREV_COL)), .before = 1)

alpha_loc_avg <- alpha_df %>%
  group_by(!!rlang::sym(LOCATION_COL)) %>%
  summarize(
    n             = dplyr::n(),
    Observed_mean = mean(Observed, na.rm = TRUE),
    Observed_sd   = sd(Observed, na.rm = TRUE),
    Shannon_mean  = mean(Shannon, na.rm = TRUE),
    Shannon_sd    = sd(Shannon, na.rm = TRUE),
    Simpson_mean  = mean(Simpson, na.rm = TRUE),
    Simpson_sd    = sd(Simpson, na.rm = TRUE),
    Chao1_mean    = mean(Chao1, na.rm = TRUE),
    Chao1_sd      = sd(Chao1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(!!LOCATION_COL := factor(.data[[LOCATION_COL]], levels = levels(sd_ps[[LOCATION_COL]])))

# Save tables
readr::write_csv(alpha_df,      file.path(outdir, "alpha_per_sample.csv"))
readr::write_csv(alpha_loc_avg, file.path(outdir, "alpha_location_averages.csv"))

plot_alpha_metric <- function(metric, ylab_txt){
  ggplot(alpha_df,
         aes(x = .data[[LOCATION_COL]],
             y = .data[[metric]],
             color = .data[[ABBREV_COL]],
             shape = .data[[ZONE_COL]])) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.6, alpha = 0.95) +
    ggrepel::geom_text_repel(aes(label = .data[[ABBREV_COL]]),
                             size = 2.4, max.overlaps = 25,
                             position = position_jitter(width = 0.12, height = 0),
                             segment.color = "grey80", alpha = 0.8, show.legend = FALSE) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 0.4) +
    stat_summary(fun.data = ggplot2::mean_se, geom = "pointrange", color = "black", linewidth = 0.3) +
    scale_color_manual(values = col_abbr, name = "Individual") +
    scale_shape_manual(values = shape_vals, name = "Zone") +
    labs(x = "Location (ordered by Zone)", y = ylab_txt,
         title = paste("Alpha diversity:", metric),
         subtitle = "Points = samples (colored by individual, shaped by Zone); black = Location mean ± SE") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

ggsave(file.path(outdir, "alpha_Shannon_by_location.png"),
       plot_alpha_metric("Shannon", "Shannon index"), width = 12, height = 6, dpi = 140)
ggsave(file.path(outdir, "alpha_Chao1_by_location.png"),
       plot_alpha_metric("Chao1", "Chao1 richness"), width = 12, height = 6, dpi = 140)
ggsave(file.path(outdir, "alpha_Observed_by_location.png"),
       plot_alpha_metric("Observed", "Observed ASVs"), width = 12, height = 6, dpi = 140)
ggsave(file.path(outdir, "alpha_Simpson_by_location.png"),
       plot_alpha_metric("Simpson", "Simpson index"), width = 12, height = 6, dpi = 140)

## =========================
## BETA DIVERSITY (Bray–Curtis) + Location centroids
## =========================
ps_rel <- transform_sample_counts(ps_use, function(x) x / sum(x))
OTU    <- as(otu_table(ps_rel), "matrix")
if (taxa_are_rows(ps_rel)) OTU <- t(OTU)  # rows = samples

d_bray <- vegan::vegdist(OTU, method = "bray")
pcoa   <- cmdscale(d_bray, k = 2, eig = TRUE)

axes <- as.data.frame(pcoa$points)
axes$SampleID <- rownames(axes)
eig_pos <- pcoa$eig[pcoa$eig > 0]
var_explained <- round(100 * (pcoa$eig[1:2] / sum(eig_pos)), 1)

axes <- axes %>%
  rename(PCoA1 = V1, PCoA2 = V2) %>%
  left_join(sd_ps %>% tibble::rownames_to_column("SampleID"), by = "SampleID") %>%
  mutate(
    !!LOCATION_COL := factor(.data[[LOCATION_COL]], levels = levels(sd_ps[[LOCATION_COL]])),
    !!ZONE_COL     := factor(.data[[ZONE_COL]],     levels = zone_levels),
    !!ABBREV_COL   := factor(.data[[ABBREV_COL]],   levels = abbr_levels)
  )

# Location centroids (average of samples per Location)
centroids <- axes %>%
  group_by(!!rlang::sym(LOCATION_COL)) %>%
  summarize(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop") %>%
  left_join(loc_by_zone %>% select(all_of(c(LOCATION_COL, ZONE_COL))), by = LOCATION_COL) %>%
  mutate(!!ZONE_COL := factor(.data[[ZONE_COL]], levels = zone_levels))

# Links: sample -> its Location centroid
links <- axes %>%
  select(SampleID, all_of(c(LOCATION_COL, "PCoA1", "PCoA2"))) %>%
  left_join(centroids %>% rename(PCoA1_c = PCoA1, PCoA2_c = PCoA2), by = LOCATION_COL)

# Plot (FIXED geom_segment: provide start AND end)
p_beta <- ggplot(axes,
                 aes(x = PCoA1, y = PCoA2,
                     color = .data[[ABBREV_COL]],
                     shape = .data[[ZONE_COL]])) +
  # sample → centroid links
  geom_segment(data = links,
               aes(x = PCoA1, y = PCoA2, xend = PCoA1_c, yend = PCoA2_c),
               inherit.aes = FALSE, linewidth = 0.3, color = "grey70") +
  # sample points
  geom_point(size = 2.6, alpha = 0.95) +
  # centroid points (bigger, outlined)
  geom_point(data = centroids,
             aes(x = PCoA1, y = PCoA2, shape = .data[[ZONE_COL]]),
             inherit.aes = FALSE, size = 4, color = "black", stroke = 1) +
  # centroid labels (need x & y here!)
  ggrepel::geom_text_repel(
    data = centroids,
    aes(x = PCoA1, y = PCoA2, label = .data[[LOCATION_COL]]),
    inherit.aes = FALSE, size = 3.1,
    box.padding = 0.5, min.segment.length = 0
  ) +
  # sample labels (inherits x & y from ggplot)
  ggrepel::geom_text_repel(
    aes(label = .data[[ABBREV_COL]]),
    size = 2.2, max.overlaps = 30, show.legend = FALSE
  ) +
  scale_color_manual(values = col_abbr, name = "Individual") +
  scale_shape_manual(values = shape_vals, name = "Zone") +
  labs(
    title = "PCoA of Bray–Curtis with Location centroids",
    subtitle = "Points = samples (colored by individual, shaped by Zone). Large outlined points = Location averages.",
    x = paste0("PCoA 1 (", var_explained[1], "%)"),
    y = paste0("PCoA 2 (", var_explained[2], "%)")
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "beta_pcoa_bray_with_location_centroids.png"),
       p_beta, width = 10.5, height = 7.5, dpi = 140)

# Save coordinates
readr::write_csv(axes,      file.path(outdir, "beta_pcoa_bray_points.csv"))
readr::write_csv(centroids, file.path(outdir, "beta_pcoa_bray_centroids_by_location.csv"))

## =========================
## Differential Abundance (ANCOMBC, not ANCOMBC2)
## =========================
# ---- ASV level ----
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ANCOMBC)
})

outdir <- "exports"; dir.create(outdir, showWarnings = FALSE)

# ---------- Choose a reference level for Zone ----------
sd <- sample_data(ps_use)
sd$Zone <- factor(sd$Zone)
# Use most frequent Zone as reference (or set manually: ref <- "YourRef")
ref <- names(sort(table(sd$Zone), decreasing = TRUE))[1]
sd$Zone <- stats::relevel(sd$Zone, ref = ref)
sample_data(ps_use) <- sd
message("ANCOM-BC reference level for Zone: ", ref)

# ---------- Helper to tidy ANCOM-BC result ----------
tidy_ancombc <- function(res, ps_obj){
  mats <- c("beta","se","W","p_val","q_val","diff_abn")
  dfs <- lapply(mats, function(m){
    x <- res[[m]]
    if (is.null(x)) return(NULL)
    df <- as.data.frame(x); df$taxa <- rownames(x)
    pivot_longer(df, -taxa, names_to = "coef", values_to = m)
  })
  out <- Reduce(function(a,b) dplyr::full_join(a,b, by=c("taxa","coef")), dfs)
  taxdf <- as.data.frame(tax_table(ps_obj)); taxdf$taxa <- rownames(taxdf)
  left_join(out, taxdf, by = "taxa")
}

# =====================================================
# ASV level
# =====================================================
res_asv <- ANCOMBC::ancombc(
  phyloseq     = ps_use,
  formula      = "Zone",
  p_adj_method = "holm",
  prv_cut      = 0.10,
  lib_cut      = 1000,
  group        = "Zone",
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  tol          = 1e-5,
  max_iter     = 100,
  conserve     = TRUE,   # small-sample safeguard; replaces the idea behind s0_perc
  alpha        = 0.05,
  global       = TRUE
)

asv_tbl <- tidy_ancombc(res_asv$res, ps_use)
asv_sig <- asv_tbl %>% filter(q_val < 0.05)

write.csv(asv_tbl, file.path(outdir, "DA_ANCOMBC_ASV_full.csv"), row.names = FALSE)
write.csv(asv_sig, file.path(outdir, "DA_ANCOMBC_ASV_sig_q<0.05.csv"), row.names = FALSE)

# Volcano per contrast at ASV level
if (!is.null(res_asv$res$beta)) {
  for (cf in colnames(res_asv$res$beta)) {
    dfc <- asv_tbl %>% filter(coef == cf) %>% mutate(sig = q_val < 0.05)
    p <- ggplot(dfc, aes(x = beta, y = -log10(p_val), color = sig)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("FALSE"="#999999","TRUE"="#d62728")) +
      labs(title = paste0("ANCOM-BC (ASV): ", cf, " vs ", ref),
           x = "log fold-change (beta)", y = "-log10(p)") +
      theme_minimal(base_size = 12) + theme(legend.position = "none")
    ggsave(file.path(outdir, paste0("DA_ANCOMBC_ASV_volcano_", cf, ".png")),
           p, width = 8, height = 6, dpi = 140)
  }
}

# =====================================================
# Genus level (aggregate first, then run ancombc)
# =====================================================
if ("Genus" %in% colnames(tax_table(ps_use))) {
  ps_genus <- tax_glom(ps_use, taxrank = "Genus", NArm = TRUE)
  # keep the same Zone reference
  sdg <- sample_data(ps_genus); sdg$Zone <- stats::relevel(factor(sdg$Zone), ref = ref)
  sample_data(ps_genus) <- sdg
  
  res_genus <- ANCOMBC::ancombc(
    phyloseq     = ps_genus,
    formula      = "Zone",
    p_adj_method = "holm",
    prv_cut      = 0.10,
    lib_cut      = 1000,
    group        = "Zone",
    struc_zero   = TRUE,
    neg_lb       = TRUE,
    tol          = 1e-5,
    max_iter     = 100,
    conserve     = TRUE,
    alpha        = 0.05,
    global       = TRUE
  )
  
  genus_tbl <- tidy_ancombc(res_genus$res, ps_genus)
  genus_sig <- genus_tbl %>% filter(q_val < 0.05)
  
  write.csv(genus_tbl, file.path(outdir, "DA_ANCOMBC_Genus_full.csv"), row.names = FALSE)
  write.csv(genus_sig, file.path(outdir, "DA_ANCOMBC_Genus_sig_q<0.05.csv"), row.names = FALSE)
  
  # Volcano per contrast at Genus level
  if (!is.null(res_genus$res$beta)) {
    for (cf in colnames(res_genus$res$beta)) {
      dfc <- genus_tbl %>% filter(coef == cf) %>% mutate(sig = q_val < 0.05)
      p <- ggplot(dfc, aes(x = beta, y = -log10(p_val), color = sig)) +
        geom_point(alpha = 0.8) +
        scale_color_manual(values = c("FALSE"="#999999","TRUE"="#d62728")) +
        labs(title = paste0("ANCOM-BC (Genus): ", cf, " vs ", ref),
             x = "log fold-change (beta)", y = "-log10(p)") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
      ggsave(file.path(outdir, paste0("DA_ANCOMBC_Genus_volcano_", cf, ".png")),
             p, width = 8, height = 6, dpi = 140)
    }
  }
}

####
####
## ===== Significance testing for beta diversity =====
## Requires: ps_use (phyloseq with otu_table + sample_data),
##           ps_tree (same samples + a rooted tree) for UniFrac
suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
})

stopifnot(inherits(ps_use, "phyloseq"))
md <- as(sample_data(ps_use), "data.frame")
stopifnot(all(c("Zone","Location") %in% names(md)))

outdir <- "exports"; dir.create(outdir, showWarnings = FALSE)

## --- Distance matrices ---
# Bray–Curtis on relative abundance
ps_rel  <- transform_sample_counts(ps_use, function(x) x / sum(x))
OTU_rel <- as(otu_table(ps_rel), "matrix"); if (taxa_are_rows(ps_rel)) OTU_rel <- t(OTU_rel)
dist_bray <- vegdist(OTU_rel, method = "bray")

# Jaccard on presence/absence (non-phylogenetic)
OTU_pa <- (as(otu_table(ps_use), "matrix") > 0) * 1
if (taxa_are_rows(ps_use)) OTU_pa <- t(OTU_pa)
dist_jaccard <- vegdist(OTU_pa, method = "jaccard", binary = TRUE)

# UniFrac (needs ps_tree with a rooted tree)
dist_unifrac_unw <- phyloseq::distance(ps_tree, method = "unifrac")   # unweighted
dist_unifrac_w   <- phyloseq::distance(ps_tree, method = "wunifrac")  # weighted

## --- Helper to run PERMANOVA + betadisper and save text reports ---
run_tests <- function(D, D_name, factor_var = c("Zone","Location"), perms = 999){
  for (f in factor_var){
    message(sprintf("Testing %s by %s ...", D_name, f))
    # PERMANOVA
    perm <- adonis2(D ~ md[[f]], data = md, permutations = perms)
    cap_perm <- capture.output(
      cat(sprintf("\n=== PERMANOVA: %s ~ %s ===\n", D_name, f)),
      print(perm)
    )
    writeLines(cap_perm, file.path(outdir, sprintf("permanova_%s_by_%s.txt", D_name, f)))
    # Homogeneity of group dispersion
    bd <- betadisper(D, md[[f]])
    cap_bd <- capture.output(
      cat(sprintf("\n=== betadisper ANOVA: %s by %s ===\n", D_name, f)),
      print(anova(bd)),
      cat("\n--- Permutation test of dispersion ---\n"),
      print(permutest(bd, permutations = perms)),
      cat("\n--- Pairwise Tukey HSD on dispersion ---\n"),
      print(TukeyHSD(bd))
    )
    writeLines(cap_bd, file.path(outdir, sprintf("betadisper_%s_by_%s.txt", D_name, f)))
  }
}

## --- Run tests and write text files into exports/ ---
run_tests(dist_bray,         "BrayCurtis", c("Zone","Location"))
run_tests(dist_jaccard,      "Jaccard",    c("Zone","Location"))
run_tests(dist_unifrac_unw,  "UniFrac_unweighted", c("Zone","Location"))
run_tests(dist_unifrac_w,    "UniFrac_weighted",   c("Zone","Location"))

message("Done. See text reports in 'exports/' for PERMANOVA tables and dispersion diagnostics.")





