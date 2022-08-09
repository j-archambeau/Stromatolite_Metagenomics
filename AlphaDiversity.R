# set directory to your working drive
setwd("~/REU/")

# upload necessary libraries
library(SRS)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(vegan)
library(rlang)
library(gt)
library(readr)
library(patchwork)
library(zCompositions)
library(mixOmics)
library(ggpubr)
library(corrplot)
library(psych)
library(ggforce)

##                         Alpha Diversity

# import biom file
physeq <- import_biom("cuatrocAssemblies.biom")

#make the phyloseq object more concise
physeq@tax_table@.Data <- substring(physeq@tax_table@.Data, 4) 
  # might have to run twice to remove unnecessary labels
colnames(physeq@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", 
                                      "Order", "Family", "Genus", 
                                      "Species")
# filter by including only bacteria
unique(physeq@tax_table@.Data[,"Phylum"])
physeq <- subset_taxa(physeq, Kingdom == "Bacteria")

#get rid of any uncharacterized OTUs down to Genus level
summary(physeq@tax_table@.Data=="")
physeq <- subset_taxa(physeq, Genus !="")
summary(physeq@tax_table@.Data=="")      

#normalize OTU counts with SRS package
otu_table.df <- as.data.frame(otu_table(physeq))
table(colSums(otu_table.df))

#Cmin = lowest number of OTU from table made above
# double check value with equation
Cmin <- min(colSums(otu_table.df))
Cmin

SRS.scaled.otu<-SRS(otu_table.df, Cmin=360)

#check that the normalization worked
p.raw.curve <- SRScurve(otu_table.df, "shannon")
p.SRS.curve <-SRScurve(SRS.scaled.otu, "shannon")
p.raw.curve                            
p.SRS.curve

#make the scaled OTU dataframe into an otu_table for phyloseq
scaled.otu_table <- otu_table(SRS.scaled.otu,taxa_are_rows=TRUE)

#estimate richness & plot richness (chose what best answers the questions
  # you're asking)
SRS.diversity<- estimate_richness(scaled.otu_table, split=TRUE, 
                                  measures = c("Shannon", "Simpson",
                                               "Chao1"))
SRS.diversity
row.names(SRS.diversity)[1] <- "Library 1"
row.names(SRS.diversity)[2] <- "Library 2"
row.names(SRS.diversity)[3] <- "Library 3"
row.names(SRS.diversity)[4] <- "Library 4"
colnames(SRS.diversity)[2] <- "Chao1 SE"
scaled.otu_table

plot_rich<-plot_richness(physeq = scaled.otu_table,
                         measures = c("Shannon", "Simpson", "Chao1"),
                         title = 
                         "Alpha Diversity Measurements on Normalized 
            Stromatolite Assemblies") 
plot_rich

tick_marks <- c("Library 1", "Library 2", "Library 3", "Library 4")

Richness <- plot_rich + scale_x_discrete(
  guide = guide_axis(angle = 90), labels = tick_marks) +
  scale_y_continuous(limits = c(0,1000), c(0,5), c(0,1))
Richness
