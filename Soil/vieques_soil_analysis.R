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
  coord_flip()+
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
#31.2
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#12.5

#plot all samples
ggplot(ko.coords, aes(MDS1, MDS2, label=SampleID, color=Soil))+
  geom_point(size=2.7)+
  #geom_text()+
  theme_bw()+
  ylab("PC2- 12.2%")+
  guides(alpha = "none")+
  xlab("PC1- 31.2%")

#beta diversity for erach site only
meta_com<-meta[which(meta$Soil == 'Composta'),]
meta_sem<-meta[-which(meta$Soil == 'Composta'),]

otu_comp<-asv.tbl[,names(asv.tbl) %in% meta_com$SampleID]
otu_sem<-asv.tbl[,names(asv.tbl) %in% meta_sem$SampleID]

#rarefy and calculate beta div
comp_pcoa<-capscale(t(otu_comp) ~1, distance = 'bray')
sem_pcoa<-capscale(t(otu_sem) ~1, distance='bray')

#pull out x/y coordinates
sem.scores<-scores(sem_pcoa)
comp.scores<-scores(comp_pcoa)

#grab only sample coordinates, write to data frame
sem.coords<-as.data.frame(sem.scores$sites)
comp.coords<-as.data.frame(comp.scores$sites)

#create sample names as a column
sem.coords$SampleID<-row.names(sem.coords)
comp.coords$SampleID<-row.names(comp.coords)

#map back meta data
sem.coords<-merge(sem.coords, meta, by.x='SampleID', by.y='SampleID', all.y=F)
comp.coords<-merge(comp.coords, meta, by.x='SampleID', by.y='SampleID', all.y=F)

ggplot(sem.coords, aes(MDS1, MDS2, label=SampleID, color=PlantSpecies))+
  geom_point(size=2.7)+
  #geom_text()+
  theme_bw()+
  ylab("PC2- 12.2%")+
  ggtitle("La Semillera")+
  guides(alpha = "none")+
  xlab("PC1- 31.2%")

ggplot(comp.coords, aes(MDS1, MDS2, label=SampleID, color=PlantSpecies))+
  geom_point(size=2.7)+
  ggtitle('Conposta')+
  #geom_text()+
  theme_bw()+
  ylab("PC2- 12.2%")+
  guides(alpha = "none")+
  xlab("PC1- 31.2%")

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

