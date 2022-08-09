# create taxonomic bar plots by phylum and genus for Anza Borrego stromatolite reads
str_metagenomes <- import_biom("cuatrocReads.biom")
str_metagenomes@tax_table@.Data <- substring(str_metagenomes@tax_table@.Data, 4) 
View(str_metagenomes@tax_table@.Data)
# might have to run line 357 twice to remove unnecessary characters

colnames(str_metagenomes@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", 
                                      "Order", "Family", "Genus", 
                                      "Species")
#add column names to subset
str_metagenomes <- subset_taxa(str_metagenomes, Genus != "")
View(str_metagenomes@otu_table@.Data)
plotnames <- c("Library 1", "Library 2",
                "Library 3", "Library 4")
View(str_metagenomes@otu_table@.Data)
colnames(str_metagenomes@otu_table@.Data)<- plotnames

#change abundance to percentages
percentages_str  = transform_sample_counts(str_metagenomes, function(x) x*100 / sum(x) )
str_glom_p <- tax_glom(percentages_str, taxrank = 'Phylum')
str_glom_g <- tax_glom(percentages_str, taxrank = 'Genus')
str_phyla <- psmelt(str_glom_p)
str_genus <- psmelt(str_glom_g)
str_phyla$Phylum[str_phyla$Abundance < 0.5] <- "Phyla < 0.5% abund."
str_genus$Genus[str_genus$Abundance < 0.5] <- "Genus < 0.5% abund."
str_phyla$Sample <- factor(str_phyla$Sample, levels=plotnames)
str_genus$Sample <- factor(str_genus$Sample, levels=plotnames)
str.phy.plot <- ggplot(data=str_phyla, aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack")
str.gen.plot <- ggplot(data=str_genus, aes(x=Sample, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack")
str.phy.plot
str.gen.plot
