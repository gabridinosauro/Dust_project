library(vegan)
library(reshape2)



####Amount of NAs, dust vs soil
asvtab<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
####add a column for sample repetitions (merge sample and repetition)
sam_data = data.frame(asvtab@sam_data)
sam_data = sam_data[,-1] #eliminate first column (useless)
sam_data$sample <- paste(sam_data$Point_type,sam_data$Repetition, sep = "_")
#####insert new asvtab into the phyloseq object
asvtab@sam_data =sample_data(sam_data) 
asvtab@sam_data$Repetition=factor(asvtab@sam_data$Repetition)
asvtab1 = data.frame(t(asvtab@otu_table))
####change the names in the graph:
sam_data$Dust_Type[sam_data$Dust_Type == "Dust_Gen"] <- "WT"
sam_data$Dust_Type[sam_data$Dust_Type == "500"] <- "Coarse SS"
sam_data$Dust_Type[sam_data$Dust_Type == "75"] <- "Medium SS"
sam_data$Dust_Type[sam_data$Dust_Type == "25"] <- "Fine SS"

sam_data$Dust_Type<- factor(sam_data$Dust_Type,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))
taxonomy = data.frame(asvtab@tax_table)
rownames(sam_data) = colnames(asvtab1)


####Average abundances
asvtab1 = decostand(asvtab1, method = "total", MARGIN = 2)

################species##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Species = taxonomy$Species
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Species +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Species~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###reorder the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = FSA::dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
rcompanion::cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_species = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + ylab("Relative abundance") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(geom = 'text', label = c("b","a","a","a","a" ), fun = max, vjust = 0) 
unknown_species




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Species == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_species_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("b","a","a","ab","a" ), fun = max, vjust = 0)
unknown_species_freq

unknown_arranged_species = ggarrange(unknown_species,unknown_species_freq, ncol = 1 )
unknown_arranged_species = annotate_figure(unknown_arranged_species, top = text_grob("Unassigned species ASVs"))

unknown_arranged_species









################Genus##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Genus = taxonomy$Genus
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Genus +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Genus~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###reorder the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Genus = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + ylab("Relative abundance") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(geom = 'text', label = c("b","a","a","a","a" ),fun = max, vjust = 0)# + ggtitle ("Bact/Arch, unassigned genus")
unknown_Genus
#saveRDS(unknown_Genus, "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4panelG.RDS")
ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4/panel123.pdf",   device = pdf)




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Genus == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Genus_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("b","a","a","ab","a" ), fun = max, vjust = 0)
unknown_Genus_freq

unknown_arranged_Genus = ggarrange(unknown_Genus,unknown_Genus_freq, ncol = 1 )
unknown_arranged_Genus = annotate_figure(unknown_arranged_Genus, top = text_grob("Unassigned Genus ASVs"))

unknown_arranged_Genus










################Family##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Family = taxonomy$Family
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Family +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Family~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###reorder the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Family = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("b","a","a","a","a" ), fun = max, vjust = 0)
unknown_Family




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Family == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Family_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency") + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("b","a","a","ab","a" ), fun = max, vjust = 0)
unknown_Family_freq

unknown_arranged_Family = ggarrange(unknown_Family,unknown_Family_freq, ncol = 1 )
unknown_arranged_Family = annotate_figure(unknown_arranged_Family, top = text_grob("Unassigned Family ASVs"))
unknown_arranged_Family




################Order##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Order = taxonomy$Order
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Order +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Order~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###reorder the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Order = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("b","a","a","a","a" ), fun = max, vjust = 0)
unknown_Order




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Order == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Order_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency") + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("b","a","a","ab","a" ), fun = max, vjust = 0)
unknown_Order_freq

unknown_arranged_Order = ggarrange(unknown_Order,unknown_Order_freq, ncol = 1 )
unknown_arranged_Order = annotate_figure(unknown_arranged_Order, top = text_grob("Unassigned Order ASVs"))

unknown_arranged_Order



################ CLass  ##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Class = taxonomy$Class
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Class +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Class~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###reClass the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Class = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("b","a","ab","ab","ab" ), fun = max, vjust = 0)
unknown_Class




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Class == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Class_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency") + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Class_freq

unknown_arranged_Class = ggarrange(unknown_Class,unknown_Class_freq, ncol = 1 )
unknown_arranged_Class = annotate_figure(unknown_arranged_Class, top = text_grob("Unassigned Class ASVs"))
unknown_arranged_Class


################ Phylum  ##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Phylum = taxonomy$Phylum
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Phylum +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Phylum~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###rePhylum the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Phylum = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Phylum




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Phylum == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Phylum_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency") + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Phylum_freq

unknown_arranged_Phylum = ggarrange(unknown_Phylum,unknown_Phylum_freq, ncol = 1 )
unknown_arranged_Phylum = annotate_figure(unknown_arranged_Phylum, top = text_grob("Unassigned Phylum ASVs"))
unknown_arranged_Phylum

####
a = ggarrange(unknown_arranged_species,unknown_arranged_Genus ,unknown_arranged_Family,unknown_arranged_Order, unknown_arranged_Class, unknown_arranged_Phylum, nrow = 1 )
a  = annotate_figure(a , top = text_grob("Bacteria/Archaea"))
a
ggsave( "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/unassigned_taxonomy.pdf",   
        device = pdf,
        width = 13, height = 6.5)




###########FUNGI############
####Amount of NAs, dust vs soil
asvtab<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
####add a column for sample repetitions (merge sample and repetition)
sam_data = data.frame(asvtab@sam_data)
sam_data = sam_data[,-1] #eliminate first column (useless)
sam_data$sample <- paste(sam_data$Point_type,sam_data$Repetition, sep = "_")
#####insert new asvtab into the phyloseq object
asvtab@sam_data =sample_data(sam_data) 
asvtab@sam_data$Repetition=factor(asvtab@sam_data$Repetition)
asvtab1 = data.frame(t(asvtab@otu_table))
####change the names in the graph:
sam_data$Dust_Type[sam_data$Dust_Type == "Dust_Gen"] <- "WT"
sam_data$Dust_Type[sam_data$Dust_Type == "500"] <- "Coarse SS"
sam_data$Dust_Type[sam_data$Dust_Type == "75"] <- "Medium SS"
sam_data$Dust_Type[sam_data$Dust_Type == "25"] <- "Fine SS"

sam_data$Dust_Type<- factor(sam_data$Dust_Type,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))
taxonomy = data.frame(asvtab@tax_table)
rownames(sam_data) = colnames(asvtab1)


####Average abundances
asvtab1 = decostand(asvtab1, method = "total", MARGIN = 2)

################species##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Species = taxonomy$Species
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Species +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Species~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###reorder the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_species = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + ylab("Relative abundance") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_species




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Species == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_species_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_species_freq

unknown_arranged_species = ggarrange(unknown_species,unknown_species_freq, ncol = 1 )
unknown_arranged_species = annotate_figure(unknown_arranged_species, top = text_grob("Unassigned species ASVs"))










################Genus##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Genus = taxonomy$Genus
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Genus +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Genus~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###reorder the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Genus = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + ylab("Relative abundance") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)#+ggtitle("Fungi, unassigned genus")

unknown_Genus
#saveRDS(unknown_Genus, "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4panelF.RDS")
#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4/panel1234.pdf",   device = pdf)




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Genus == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Genus_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0) 
unknown_Genus_freq

unknown_arranged_Genus = ggarrange(unknown_Genus,unknown_Genus_freq, ncol = 1 )
unknown_arranged_Genus = annotate_figure(unknown_arranged_Genus, top = text_grob("Unassigned Genus ASVs"))

unknown_arranged_Genus











################Family##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Family = taxonomy$Family
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Family +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Family~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###reorder the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Family = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter" ) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Family




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Family == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Family_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency") + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Family_freq

unknown_arranged_Family = ggarrange(unknown_Family,unknown_Family_freq, ncol = 1 )
unknown_arranged_Family = annotate_figure(unknown_arranged_Family, top = text_grob("Unassigned Family ASVs"))





################Order##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Order = taxonomy$Order
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Order +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Order~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###reorder the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Order = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Order




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Order == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Order_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency") + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Order_freq

unknown_arranged_Order = ggarrange(unknown_Order,unknown_Order_freq, ncol = 1 )
unknown_arranged_Order = annotate_figure(unknown_arranged_Order, top = text_grob("Unassigned Order ASVs"))





################ CLass  ##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Class = taxonomy$Class
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Class +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Class~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###reClass the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Class = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Class




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Class == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Class_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency") + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Class_freq

unknown_arranged_Class = ggarrange(unknown_Class,unknown_Class_freq, ncol = 1 )
unknown_arranged_Class = annotate_figure(unknown_arranged_Class, top = text_grob("Unassigned Class ASVs"))



################ Phylum  ##############
####Now aggregate/cast
asvtab2 = asvtab1
asvtab2$Phylum = taxonomy$Phylum
asvtab2[is.na(asvtab2)] <- 'unassigned'
asvtab2_melt = melt(asvtab2)
#melted_func_tab$guild = guilds[match(melted_func_tab$variable, row.names(guilds)), 1]
aggregated <- aggregate(value~ Phylum +variable,asvtab2_melt,sum)
casted = dcast(aggregated,  Phylum~ variable)



###reoroganize the rows###
rownames(casted) = casted[,1]
casted = casted[,-1]
casted =  data.table::setcolorder(casted, as.character(rownames(sam_data))) ###rePhylum the columns
casted = data.frame(t(casted[rowSums(casted[])>0,] ))###drop rows with all

##plot the unknown###
unknown_plot = cbind.data.frame(`Relative abundance` = casted$unassigned , dust_fraction = sam_data$Dust_Type)


unknown_plot$dust_fraction<- factor(unknown_plot$dust_fraction,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))





kruskal.test(unknown_plot$`Relative abundance` ~ unknown_plot$dust_fraction)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.0006251
DT = dunnTest(`Relative abundance` ~ dust_fraction,
              data=unknown_plot,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Phylum = ggplot(unknown_plot, aes(x=dust_fraction, y=`Relative abundance`)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Phylum




###relative_frequency
asvtab3 = asvtab2[which(asvtab2$Phylum == "unassigned"),]
sums = colSums(asvtab3 != 0)
total = colSums(asvtab2 != 0)

relative_frequencies = (sums/total)[1:60]
relative_frequencies = data.frame(`Relative frequency`= relative_frequencies, Sample_type = sam_data$Dust_Type )



kruskal.test(relative_frequencies$Relative.frequency ~ relative_frequencies$Sample_type)
#Kruskal-Wallis chi-squared = 19.505, df = 4, p-value = 0.00115

DT = dunnTest(Relative.frequency~ Sample_type,
              data=relative_frequencies,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



unknown_Phylum_freq = ggplot(relative_frequencies, aes(x=Sample_type, y=Relative.frequency)) +  geom_point(position = "jitter" ) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") + xlab(NULL) +ylab("Relative frequency") + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0)
unknown_Phylum_freq

unknown_arranged_Phylum = ggarrange(unknown_Phylum,unknown_Phylum_freq, ncol = 1 )
unknown_arranged_Phylum = annotate_figure(unknown_arranged_Phylum, top = text_grob("Unassigned Phylum ASVs"))

####
a = ggarrange(unknown_arranged_species,unknown_arranged_Genus ,unknown_arranged_Family,unknown_arranged_Order, unknown_arranged_Class, unknown_arranged_Phylum, nrow = 1 )
a  = annotate_figure(a , top = text_grob("Fungi"))
a
ggsave( "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/unassigned_taxonomy_fungi.pdf",   
        device = pdf,
        width = 13, height = 6.5)
