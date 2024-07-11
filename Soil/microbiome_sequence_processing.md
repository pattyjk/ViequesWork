## Processing of microbiome samples

```
source activate qiime2-2023.5

#commented out part is what I did to import/demux in QIIME2
#load raw FASTQ reads into QIIME
#qiime tools import --type EMPSingleEndSequences --input-path ./data --output-path adan_seqs.qza

#demultiplex reads
#qiime demux emp-single \
#  --i-seqs adan_seqs.qza \
 --m-barcodes-file adan_soil_map_run2.txt \
 --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences adan_demux.qza \
  --o-error-correction-details  adan-details.qza \
  --p-no-golay-error-correction 
  
  #import both runs for analysis
 qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest.txt \
  --output-path vieques-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 
  
#quality filer
qiime quality-filter q-score \
--i-demux  adan_demux.qza \
--o-filtered-sequences  adan_demux-filtered.qza \
--o-filter-stats  adandemux-filter-stats.qza \
--q-score 25

 
 #export filter stats
  qiime tools export --input-path  adan_demux-filter-stats.qza --output-path filt_stats
 
  #call ASVs with deblur
  qiime deblur denoise-16S \
  --i-demultiplexed-seqs  adan_demux-filtered.qza \
  --p-trim-length 120 \
  --p-jobs-to-start 24 \
  --o-representative-sequences  adan_rep-seqs-deblur.qza \
  --o-table  adan_table-deblur.qza \
   --o-stats  adandeblur-stats.qza
 
 #make phylogenetic tree with fasttree
 qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences adan_rep-seqs-deblur.qza \
  --output-dir phylogeny-align-to-tree-mafft-fasttree \
  --p-n-threads 24
  
 #export deblur stats
 qiime tools export --input-path adan_deblur-stats.qza --output-path deblur_stats
  
 #export rep seqs
 qiime tools export --input-path adan_rep-seqs-deblur.qza --output-path rep_seqs
  
#export tree as NWK format
qiime tools export --input-path phylogeny-align-to-tree-mafft-fasttree/tree.qza --output-path tree
 
 #export OTU table to biom then to text file
 qiime tools export --input-path adan_table-deblur.qza --output-path asv_table
 biom convert -i asv_table/feature-table.biom --to-tsv -o asv_table.txt
 
 
 #assign taxonomy with sklearn and silva database
qiime feature-classifier classify-sklearn   --i-classifier silva-138-99-515-806-nb-classifier.qza   --i-reads  adan_rep-seqs-deblur.qza   --o-classification adan_taxonomy.qza 

 #export taxonomy file
 qiime tools export --input-path taxonomy.qzv --output-path taxonomy
 
#make stacked bar visualizations
qiime taxa barplot --o-visualization taxa_plot  --m-metadata-file adan_soil_map_run2.txt  --i-taxonomy adan_taxonomy.qza --i-table adan_table-deblur.qza
```
