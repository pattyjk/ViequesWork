#load packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

#read in metadata and asv table
set.seed(515)
asv_table <- read.delim("~/GitHub/ViequesWork/Soil/asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/GitHub/ViequesWork/Soil/PR_soil_metadata.txt", header=T)

#look at sequencing depth
min(colSums(asv_table))
#3621 is the lowest good depth

#rarefy data 
nut_rare<-rrarefy(t(asv_table), sample=3621)

#calculate PCoA based on BC similarity
ko_pcoa<-capscale(nut_rare  ~ 1, distance='bray')

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
#27.4
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#7.8

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, color=PlantSpecies, shape=Soil, label=SampleID))+
  geom_point(size=3.5)+
  #geom_text()+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 27.4%")+
  ylab("PC2- 7.8%")


#calculate PERMANOVA with adonis2
s16.dis<-as.data.frame(as.matrix(vegdist(t(asv_table), method='bray')))
s16.dis2<-s16.dis
s16.dis2$SampleID<-row.names(s16.dis2)

#select only columns of interest for the metadata
meta2<-meta[,c(1,5,8)]

#add metadata 
s16.dis3<-merge(s16.dis2, meta2, by='SampleID')

#remove sample id column
s16.dis4<-s16.dis3[,-1]

#calculate adonis
set.seed(4)
adonis2(s16.dis4[,-c(52,53)] ~ PlantSpecies + Soil, data=s16.dis4, permutations = 10000)

#             Df SumOfSqs      R2       F    Pr(>F)    
#PlantSpecies  8  0.51920 0.77545 19.7213 9.999e-05 ***
#Soil          1  0.01542 0.02303  4.6851    0.051   
#Residual     41  0.13492 0.20152                      
#Total        50  0.66954 1.00000                      

###Calculate alpha diversity

#load in asv table and metadata
asv_table <- read.delim("~/GitHub/ViequesWork/Soil/asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/GitHub/ViequesWork/Soil/PR_soil_metadata.txt", header=T)

#CALCULATE RICHNESS & add metadata & statistics
larv.alph<-as.data.frame(specnumber(rrarefy(t(asv_table), sample=3621)))
larv.alph$SampleID<-row.names(larv.alph)
larv.alph<-merge(larv.alph, meta, by='SampleID')


larv.alpha2<-as.data.frame(vegan::diversity(rrarefy(t(asv_table), sample=3621), index = 'shannon'))
names(larv.alpha2)<-"Shannon"
larv.alph<-cbind(larv.alph, larv.alpha2)


#plot richness
ggplot(larv.alph, aes(PlantSpecies, `specnumber(rrarefy(t(asv_table), sample = 3621))`, fill=PlantSpecies))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  guides(fill="none")+
  ylab("sOTU Richness")

ggplot(larv.alph, aes(PlantSpecies, Shannon, fill=PlantSpecies))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("Shannon Diversity")

pairwise.t.test(larv.alph$Shannon, larv.alph$Soil)
#             Composta
#La Semillera 0.092

pairwise.t.test(larv.alph$Shannon, larv.alph$PlantSpecies)
