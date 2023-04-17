# set directory to your working drive
setwd("~/Desktop/REU/merged_metagenome")

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

### Color palettes ###
cbPalette1 <- c("#D55E00", "#56B4E9", "#9944A9", "#6511B4","#999999", 
                         "#0072B2", "#000000", "#661100", "#CC79A7","#DDCC77", "#117733", "#332288")
                         
cbPalette2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                         "#661100", "#6511B4", "#000000")
                         
cbPalette3 <- c("#56B4E9", "#9944A9", "#6511B4","#999999", 
                         "#0072B2", "#000000", "#661100", "#CC79A7","#DDCC77", "#117733", "#332288")
                         
cbPalette4 <- c("#E69F00", "#009E73", "#661100", "#6511B4", "#000000")

cbPalette5 <- c("#D55E00", "#56B4E9", "#6511B4", "#117733")

cbPalette6 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                         "#661100", "#6511B4", "#000000")

### Reads -- Comparing Metagenomes ###

## Alpha diversity ##
reads <- import_biom("full_scaffold_beta.biom")
reads@tax_table@.Data <- substring(reads@tax_table@.Data, 4)
# may need to run line 42 twice -- check file for excess "_"
colnames(reads@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", 
                                    "Order", "Family", "Genus", 
                                    "Species")
summary(reads@tax_table@.Data== "")
reads_beta <- subset_taxa(reads, Genus != "") # remove blanks
summary(reads_beta@tax_table@.Data== "")
otu_reads.df <- as.data.frame(otu_table(reads_beta))

#Cmin = lowest number of OTU from table made above
# double check value with equation
Cmin_reads <- min(colSums(otu_reads.df))
Cmin_reads
SRS.scaled.otu.reads<-SRS(otu_reads.df, Cmin=3861)
#check that normalization worked
p.raw.curve.reads <- SRScurve(otu_reads.df, "shannon")
p.SRS.curve.reads <-SRScurve(SRS.scaled.otu.reads, "shannon")
p.raw.curve.reads                       
p.SRS.curve.reads
#make the scaled OTU dataframe into an otu_table for phyloseq
scaled.otu.reads_table <- otu_table(SRS.scaled.otu.reads,taxa_are_rows=TRUE)
#estimate richness & plot richness with Shannon, Simpson, and Chao1
SRS.diversity.reads<- estimate_richness(scaled.otu.reads_table, split=TRUE, 
                                        measures = c("Shannon", "Simpson",
                                                     "Chao1"))
colnames(SRS.diversity.reads)[2] <- "Chao1 SE"
SRS.diversity.full.reads <- SRS.diversity.reads[-c(2:5),]
SRS.diversity.full.reads

## Taxonomic bar plots ##
percentages_reads = transform_sample_counts(reads_beta, function(x) x*100/sum(x)) # relative abundances calculated by taking percentages
plotnames <- c("Anza Borrego","Bahamas", "CSharkBay", "NSharkBay", "SSharkBay")
colnames(percentages_reads@otu_table@.Data)<- plotnames
sample_names(percentages_reads) <- factor(sample_names(percentages_reads), levels = list(full = "Anza Borrego", Bahamas_megahit = "Bahamas",
                                                                                         CSharkBay_megahit = "CSharkBay", NSharkBay_megahit = "NSharkBay",
                                                                                         SSharkBay_megahit = "SSharkBay"))
sample_names(percentages_reads)
sample_data(percentages_reads) # Make sure labels match

glom_p_r <- tax_glom(percentages_reads, taxrank = 'Phylum')
glom_g_r <- tax_glom(percentages_reads, taxrank = 'Genus')
reads.phyla <- psmelt(glom_p_r)
reads.genus <- psmelt(glom_g_r)

reads.phyla$Phylum[reads.phyla$Phylum == "Firmicutes_A"] <- "Firmicutes" #combine phyla

reads.genus$Genus[reads.genus$Genus == "PCC-7108"] <- "Anabaena"
reads.genus$Genus[reads.genus$Genus == "Anabaena_A"] <- "Anabaena"
reads.genus$Genus[reads.genus$Genus == "Nostoc_B"] <- "Nostoc"
reads.genus$Genus[reads.genus$Genus == "Nostoc_C"] <- "Nostoc" #combine genera

reads.phyla$Phylum[reads.phyla$Abundance < 2] <- "Phyla < 2% abund."
reads.genus$Genus[reads.genus$Abundance < 2] <- "Genera < 2% abund." #set limit for percent abundance
reads.phyla$Sample <- factor(reads.phyla$Sample, levels=plotnames)
reads.genus$Sample <- factor(reads.genus$Sample, levels=plotnames)

reads.phy.plot <- ggplot(data=reads.phyla, aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values=cbPalette1)
reads.phy.plot

reads.gen.plot <- ggplot(data=reads.genus, aes(x=Sample, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values=cbPalette2)
reads.gen.plot

### Assemblies -- Comparing Metagenomes ###

## Alpha diversity ##
# import biom file
merged_metagenomes <- import_biom("full_beta.biom")
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
# may need to run line 114 twice -- check file for excess "_"
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", 
                                                 "Order", "Family", "Genus", 
                                                 "Species")
summary(merged_metagenomes@tax_table@.Data== "")
merged_metagenomes_beta <- subset_taxa(merged_metagenomes, Genus != "") # remove blanks
summary(merged_metagenomes_beta@tax_table@.Data== "")
# see that the blanks have disappeared, counts are lower
head(merged_metagenomes_beta@otu_table@.Data)
# convert otus into a data frame
otu_table.df <- as.data.frame(otu_table(merged_metagenomes_beta))
#Cmin = lowest number of OTU from table made above
# double check value with equation
Cmin <- min(colSums(otu_table.df)) #do we need colSums? all in one column
Cmin
View(otu_table.df)
SRS.scaled.otu<-SRS(otu_table.df, Cmin=3861)
#check that the normalization worked
p.raw.curve <- SRScurve(otu_table.df, "shannon")
p.SRS.curve <-SRScurve(SRS.scaled.otu, "shannon")
p.raw.curve                            
p.SRS.curve
#make the scaled OTU dataframe into an otu_table for phyloseq
scaled.otu_table <- otu_table(SRS.scaled.otu,taxa_are_rows=TRUE)
#estimate richness & plot richness using Shannon, Simpson, and Chao1
SRS.diversity<- estimate_richness(scaled.otu_table, split=TRUE, 
                                  measures = c("Shannon", "Simpson",
                                               "Chao1"))
colnames(SRS.diversity)[2] <- "Chao1 SE"
SRS.diversity.full <- SRS.diversity[-c(2:5),]
SRS.diversity.full

## Taxonomic bar plots ##
#convert abundances into percentages
percentages = transform_sample_counts(merged_metagenomes_beta, function(x) x*100/sum(x))
# relative abundances calculated by taking percentages
head(percentages@otu_table@.Data)
percentages_otu <- percentages@otu_table@.Data
t(percentages_otu) #pivot table

#change sample names
plotnames <- c("Anza Borrego","Bahamas", "CSharkBay", "NSharkBay", "SSharkBay")
colnames(percentages@otu_table@.Data) <- plotnames
sample_names(percentages) <- factor(sample_names(percentages), levels = list(full = "Anza Borrego", Bahamas_megahit = "Bahamas",
                                                                             CSharkBay_megahit = "CSharkBay", NSharkBay_megahit = "NSharkBay",
                                                                             SSharkBay_megahit = "SSharkBay"))
sample_names(percentages)
sample_data(percentages) # make sure that they match

glom_p <- tax_glom(percentages, taxrank = 'Phylum')
glom_g <- tax_glom(percentages, taxrank = 'Genus')
phyla <- psmelt(glom_p)
genus <- psmelt(glom_g)

phyla$Phylum[phyla$Phylum == "Firmicutes_A"] <- "Firmicutes" #combine phyla

genus$Genus[genus$Genus == "PCC-7108"] <- "Anabaena"
genus$Genus[genus$Genus == "Anabaena_A"] <- "Anabaena"
genus$Genus[genus$Genus == "Nostoc_B"] <- "Nostoc"
genus$Genus[genus$Genus == "Nostoc_C"] <- "Nostoc" #combine genera

phyla$Phylum[phyla$Abundance < 2] <- "Phyla < 2% abund."
genus$Genus[genus$Abundance < 2] <- "Genera < 2% abund."
phyla$Sample <- factor(phyla$Sample, levels=plotnames)
genus$Sample <- factor(genus$Sample, levels=plotnames)

phy.plot <- ggplot(data=phyla, aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values=cbPalette3)
phy.plot

gen.plot <- ggplot(data=genus, aes(x=Sample, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values=cbPalette4)
gen.plot

## Comparing Reads and Assemblies on Anza Borrego ##
# By phylum
phyla.az <- subset(phyla, Sample == "Anza Borrego")
phyla.az$Source <- "Assembly"
reads.phyla.az <- subset(reads.phyla, Sample == "Anza Borrego")
reads.phyla.az$Source <- "Reads"
all.phyla.az <- rbind(phyla.az, reads.phyla.az)
phy.plot.az <- ggplot(data=all.phyla.az, aes(x=Source, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values=cbPalette5)
phy.plot.az
# By genera
genus.az <- subset(genus, Sample == "Anza Borrego")
genus.az$Source <- "Assembly"
reads.genus.az <- subset(reads.genus, Sample == "Anza Borrego")
reads.genus.az$Source <- "Reads"
all.genus.az <- rbind(genus.az, reads.genus.az)
gen.plot.az <- ggplot(data=all.genus.az, aes(x=Source, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values=cbPalette6)
gen.plot.az

### Generating Tables ###
## Cyanobacterial genera table ##
full_cyano <- filter(full_genus, Phylum == "Cyanobacteria" & Sample == "Anza Borrego")
full_cyano$Genus[full_cyano$Genus == "Tolypothrix_B"] <- "Tolypothrix"
full_cyano$Genus[full_cyano$Genus == "Tolypothrix_A"] <- "Tolypothrix"
full_cyano$Genus[full_cyano$Genus == "Tolypothrix_C"] <- "Tolypothrix" #combine genera
full_cyano[order(full_cyano$Genus),]
full_cyano_red <- full_cyano[-c(1,4:9)] #only sample, abundance, and genus
STR_Cyanobacteria_Table <- gt(data = full_cyano_red[order(full_cyano_red$Genus),]) %>%
  tab_header(title = "Assembled Stromatolite Cyanobacteria Info")
gtsave(STR_Cyanobacteria_Table, file = "STR_Cyanobacteria_Table.docx")

## Saving raw data ##
write.table(reads.phyla, "Compared Read Metagenomes Phylum Raw Data", quote =F, row.names = F, sep = "\t") ## Raw Read Compared Metagenomes phyla
write.table(reads.genus, "Compared Read Metagenomes Genera Raw Data", quote =F, row.names = F, sep = "\t") ## Raw Read Compared Metagenomes genera

write.table(phyla, "Compared Metagenomes Phylum Raw Data", quote =F, row.names = F, sep = "\t") # Assembled Compared Metagenomes phyla
write.table(genus, "Compared Metagenomes Genera Raw Data", quote =F, row.names = F, sep = "\t") # Assembled Compared Metagenomes genera

write.table(all.phyla.az, "Anza Borrego Phylum Raw Data", quote =F, row.names = F, sep = "\t") # Anza Borrego phyla
write.table(all.genus.az, "Anza Borrego Genera Raw Data", quote =F, row.names = F, sep = "\t") # Anza Borrego genera
