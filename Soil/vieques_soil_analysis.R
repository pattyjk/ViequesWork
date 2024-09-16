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
#29.9
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#12.2



ggplot(ko.coords, aes(MDS1, MDS2, color=PlantSpecies, label=SampleID, shape=Soil))+
  #geom_point(size=2.7)+
  geom_text()+
  theme_bw()+
  ylab("PC2- 12.1%")+
  guides(alpha = "none")+
  xlab("PC1- 30.0%")

#plot PCoA
ggplot(ko.coords[which(ko.coords$Soil == 'Composta'),], aes(MDS1, MDS2, color=PlantSpecies, label=SampleID))+
 # geom_point(size=2.7)+
  geom_text()+
  theme_bw()+
  ylab("PC2- 12.1%")+
  guides(alpha = "none")+
  xlab("PC1- 30.0%")+
  ylim(c(-0.5,.25))+
  xlim(c(-.75,0))

ggplot(ko.coords[-which(ko.coords$Soil == 'Composta'),], aes(MDS1, MDS2, color=PlantSpecies, label=SampleID))+
  # geom_point(size=2.7)+
  geom_text()+
  theme_bw()+
  ylab("PC2- 12.1%")+
  guides(alpha = "none")+
  xlab("PC1- 30.0%")+
  xlim(c(1,1.5))+
  ylim(c(-.75,-1.25))

#adonis time, La Sem only
meta_composta<-meta[-which(meta$Soil == 'Composta'),]
meta_composta<-meta_composta[-c(1:2),]
otu_rare2<-otu_rare[row.names(otu_rare) %in% meta_composta$SampleID,]

larv.dis<-vegdist(otu_rare2, method='bray')
larv.dis<-as.data.frame(as.matrix(larv.dis))
larv.dis$SampleID<-row.names(larv.dis)
larv.dis<-merge(larv.dis, meta_composta[,c(1,5)], by= 'SampleID', all.y=F)
larv.dis<-larv.dis[,-1]
larv.ad<-adonis(larv.dis[,-12] ~ larv.dis$PlantSpecies, permutations = 10000)
larv.ad$aov.tab

#adonis time, Composta only
meta_composta<-meta[which(meta$Soil == 'Composta'),]
otu_rare2<-otu_rare[row.names(otu_rare) %in% meta_composta$SampleID,]

larv.dis<-vegdist(otu_rare2, method='bray')
larv.dis<-as.data.frame(as.matrix(larv.dis))
larv.dis$SampleID<-row.names(larv.dis)
larv.dis<-merge(larv.dis, meta_composta[,c(1,5)], by= 'SampleID', all.y=F)
larv.dis<-larv.dis[,-1]
larv.permanova<-adonis(larv.dis[,-15] ~ larv.dis$PlantSpecies, permutations = 10000)
larv.permanova$aov.tab

