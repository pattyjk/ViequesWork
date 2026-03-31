#ITS DADA2 Analysis from FASTQ Files (R Script)

install.packages("Rcpp")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
#load libraries 
library(Rcpp)
library(dada2)
library(dada2); packageVersion("dada2")
library(ShortRead)
library(Biostrings)

#path

# Set your working directory to where your FASTQ files are
path <- "C:/Users/valer/OneDrive/Desktop/Amplicon_P.R_Sediments/ViequesWork-main/Sediment/ITS_fastq/"# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="L001_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="L001_R2_001.fastq", full.names = TRUE))

# Extract sample names from the forward read filenames
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# Print sample names to check
print(sample.names)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Function to check if forward and reverse reads are correctly paired
check_paired_files <- function(forward_files, reverse_files) {
  if (length(forward_files) != length(reverse_files)) {
    stop("The number of forward and reverse files do not match!")
  }
  
  forward_base <- sub("_R1_001.fastq.gz", "", basename(forward_files))
  reverse_base <- sub("_R2_001.fastq.gz", "", basename(reverse_files))
  
  if (!all(forward_base == reverse_base)) {
    stop("Forward and reverse files are not correctly paired!")
  }
  
  return(TRUE)
}

# Check if the files are correctly paired
check_paired_files(fnFs, fnRs)

# Define the paths for filtered output
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="L001_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="L001_R2_001.fastq", full.names = TRUE))

# Extract sample names from the forward read filenames
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Ensure the output directory exists
if (!dir.exists(filt_path)) {
  dir.create(filt_path)
}

# Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     trimLeft=30, maxN=0, maxEE=c(2,2), 
                     rm.phix=TRUE, multithread=TRUE)


# View the output
print(out)


# Learning error rates
errF <- learnErrors(filtFs, multithread=TRUE)
#61823109 total bases in 480405 reads from 29 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#62104360 total bases in 480405 reads from 29 samples will be used for learning the error rates

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#dereplication 
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#denoising

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


#merge pair reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)


#sequence table /remove chimeras 
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

#write to file
asv.table<-as.data.frame(seqtab.nochim)
write.table(asv.table, 'asv.table.txt', row.names=T, quote=F, sep='\t')

#load ASV table 
asv.table<-read.delim("C:/Users/valer/OneDrive/Desktop/Amplicon_P.R_Sediments/ViequesWork-main/Sediment/ITS_fastq/ITS_ampl_sed". header=T, row.names=1)

#load libraries 
library(vegan)

#what's lowest seq depth
min(rowSums(asv.table))
#1015

#rarefy data
asv.rare<-rrarefy(asv.table, sample=1015)

#calculate richness
alpha<-as.data.frame(specnumber(asv.rare))
alpha$SampleID<-row.names(alpha)

#calculate shannon diversity
alpha2<-as.data.frame(diversity(asv.rare, index = 'shannon' ))
alpha2$SampleID<-row.names(alpha2)

#bind shannon and richness
alpha3<-merge(alpha2, alpha, by='SampleID')
names(alpha3)<-c('SampleID', 'Shannon', 'Richness')
write.table(alpha3, 'alphadiversity.txt', sep='\t', quote=F)

#look at beta diversity
fung.pcoa<-capscale(asv.rare~1, distance = 'bray')
fung.scores<-scores(fung.pcoa)
fung.scores<-as.data.frame(fung.scores$sites)
fung.scores$SampleID<-row.names(fung.scores)
write.table(fung.scores, 'fungal_pcoa.txt', sep='\t', quote=F)

#assing tax for ITS
unite_ref <- "sh_general_release_dynamic_s_all_04.02.2020.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite_ref, multithread = TRUE, tryRC = TRUE)

#SPECIES 
taxa <- addSpecies(taxa, unite_ref)

#Output 
write.csv(seqtab.nochim, "ITS_seqtab.csv")
write.csv(taxa, "ITS_taxonomy.csv")


