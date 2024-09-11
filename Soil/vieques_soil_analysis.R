##Vieques soil code
library(ggplot2)
library(vegan)
library(ggpubr)

set.seed(515)
#load in asv table
asv.tbl<-read.delim('~/Documents/GitHub/ViequesWork/Soil/asv_table.txt', row.names = 1, header=T)

#read in simplified map
meta<-read.delim('~/Documents/GitHub/ViequesWork/Soil/PR_soil_metadata.txt', header=T)

#sequencing depth
min(colSums(asv.tbl))
#rarefy to 4100

#rarefy data to 4100, nice like slice :)
otu_rare<-rrarefy(t(asv.tbl), sample=4100)

#CALCULATE RICHNESS & add metadata & statistics
vieq.alph<-as.data.frame(diversity(otu_rare, index = 'shannon'))
vieq.alph$SampleID<-row.names(vieq.alph)
vieq.alph<-merge(vieq.alph, meta, by='SampleID')

ggplot(vieq.alph, aes(PlantSpecies, `diversity(otu_rare, index = "shannon")`, fill=Soil))+
  geom_boxplot()+
  theme_bw()+
  ylab('Shannon Diversity')+
  xlab("")

##########################################################################
############################Beta diversity################################
##########################################################################

#calculate beta diversity
ko_pcoa<-capscale(otu_rare  ~ 1, distance='bray')

#pull out x/y coordinates
ko.scores<-scores(ko_pcoa)

#grab only sample coordinates, write to data frame
ko.coords<-as.data.frame(ko.scores$sites)

#create sample names as a column
ko.coords$SampleID<-row.names(ko.coords)

#map back meta data
ko.coords<-merge(ko.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(ko_pcoa$CA$eig[1]/sum(ko_pcoa$CA$eig), 3)
#30
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#12.1

#color=Time,

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, shape=PlantSpecies, color=Soil, label=SampleID))+
 # geom_point(size=2.7)+
  geom_text()+
  theme_bw()+
  ylab("PC2- 12.1%")+
  guides(alpha = "none")+
  xlab("PC1- 30.0%")
