###Find most common genera among generalists and Specialists
library(phyloseq)
library(dplyr)
library(reshape2)
library(ggplot2)
library(vegan)
library(FSA)
library(rcompanion)
library(ggpubr)
library(gridExtra)










#####test different thresholds for community composition#######
####High abundance and High occupancy idea######
geodata  = readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_data.RDS")
ASVtab_dust<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
ASVtab_dust=prune_samples(sample_data(ASVtab_dust)$point!=26, ASVtab_dust)
ASVtab_dust <- prune_taxa(taxa_sums(ASVtab_dust) > 0, ASVtab_dust)
sam_data_dust = data.frame(ASVtab_dust@sam_data) #extract dust data from phyloseq object
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "Dust_Gen"] <- "WT"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "500"] <- "Coarse SS"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "75"] <- "Medium SS"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "25"] <- "Fine SS"
sam_data_dust$Dust_Type<- factor(sam_data_dust$Dust_Type,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))











asvtab_dust = data.frame(ASVtab_dust@otu_table) #extract ASVt table not transformed
####set a threshold for the first analysis
threshold = 1














###select the ASVs based on high occupancy
selection = round((nrow(geodata)/100)*threshold) ###define the number of ASVs to pick
high_occupancy_ASVs = geodata[order(-geodata$OC_points_percent),] ### order the table
high_occupancy_ASVs = high_occupancy_ASVs[1:selection,] ### extract those with the highest occupancy
min = min(high_occupancy_ASVs$OC_points_percent) ####correct for those having the same occupancy
high_occupancy_ASVs = geodata[which(geodata$OC_points_percent >= min),] ##finally select them











Log2foldchanges = readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/log2foldchanges_bac_WT.RDS")
Log2foldchanges_occupancy  = Log2foldchanges[rownames(Log2foldchanges) %in% rownames(high_occupancy_ASVs), ]
Log2foldchanges_occupancy$Gen_spec = "Generalists"
Log2foldchanges_occupancy$abundance = high_occupancy_ASVs[match(rownames(high_occupancy_ASVs),rownames(Log2foldchanges_occupancy)),6 ]
Log2foldchanges_occupancy = head(Log2foldchanges_occupancy[order(Log2foldchanges_occupancy$abundance, decreasing=TRUE), ], 20)











###select the ASVs based on high abundance
`%nin%` = Negate(`%in%`)
geodata_specialists = geodata[rownames(geodata) %nin% rownames(high_occupancy_ASVs), ]
selection = round((nrow(geodata_specialists)/100)*threshold) ###define the number of ASVs to pick
high_abundance_ASVs = geodata_specialists[order(-geodata_specialists$Ab_bypoint),] ### order the table
high_abundance_ASVs = high_abundance_ASVs[1:selection,] ### extract those with the highest occupancy
min_ab = min(high_abundance_ASVs$Ab_bypoint) ####correct for those having the same occupancy
high_abundance_ASVs = geodata_specialists[which(geodata_specialists$Ab_bypoint >= min_ab),] ##finally select them



Log2foldchanges_abundance  = Log2foldchanges[rownames(Log2foldchanges) %in% rownames(high_abundance_ASVs), ] #select those found in high abundance
#now we need to assess how much are they present in dust samples
dust_samples_tot = colSums(asvtab_dust)
Log2foldchanges_abundance$abundance_dust = dust_samples_tot[match(rownames(Log2foldchanges_abundance),names(dust_samples_tot))] ####associate it to the ASVs
#Log2foldchanges_abundance = head(Log2foldchanges_abundance[order(Log2foldchanges_abundance$abundance_dust, decreasing=TRUE), ], 10) #select the 10 most abundant
Log2foldchanges_abundance$Gen_spec = "Specialists"
Log2foldchanges_abundance$abundance = Log2foldchanges_abundance[match(rownames(Log2foldchanges_abundance),rownames(high_abundance_ASVs)),6 ]
Log2foldchanges_abundance = head(Log2foldchanges_abundance[order(Log2foldchanges_abundance$abundance, decreasing=TRUE), ], 20)
Log2foldchanges_abundance = Log2foldchanges_abundance[,-13]






final_log2foldchange = rbind(Log2foldchanges_occupancy, Log2foldchanges_abundance)
final_log2foldchange$dust_soil = ifelse(final_log2foldchange$log2FoldChange > 0, "Dust" , "Soil")
final_log2foldchange$name = paste(rownames(final_log2foldchange), final_log2foldchange$Genus, sep=" - ")
species = ifelse(is.na(final_log2foldchange$Species) , "sp." , final_log2foldchange$Species)
final_log2foldchange$name = paste(final_log2foldchange$name, species, sep=" ")














fine_bact = ggplot(final_log2foldchange, 
       aes(log2FoldChange, y = name, fill = dust_soil)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV")+ggtitle("Bact/Arch") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_flip()+
  scale_fill_manual(name = "Sample", labels = c("Dust", "Soil"), values = c("4DBBD5B2","#996600")) + ylab(NULL)+
  facet_wrap(~ Gen_spec, scales = "free_x")
fine_bact


fine_bact

#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/3rd_draft/plots/new_figures/bag_gen_spec.pdf", height = 4, width = 8)

#dev.off()
















#Associate genera to ASVs
taxa_dust = data.frame(ASVtab_dust@tax_table) #extract dust data from phyloseq object
high_occupancy_ASVs$genus = taxa_dust[match(rownames(high_occupancy_ASVs), rownames(taxa_dust)), ]
high_occupancy_ASVs$genus[is.na(high_occupancy_ASVs$genus)] <- 'unknown'
#count Nuber of unknown 
a = table(high_occupancy_ASVs$genus$Genus)
 a[names(a) == "unknown"] / nrow(high_occupancy_ASVs)

high_abundance_ASVs$genus = taxa_dust[match(rownames(high_abundance_ASVs), rownames(taxa_dust)), ]
high_abundance_ASVs$genus[is.na(high_abundance_ASVs$genus)] <- 'unknown'
aggregated_ab = aggregate(Ab_total ~ genus$Genus, high_abundance_ASVs, FUN = sum)
b = table(high_abundance_ASVs$genus$Genus)
b[names(b) == "unknown"] / nrow(high_abundance_ASVs)
sequences  = readRDS("/Users/gabri/OneDrive - University of Arizona/dust_project/data/soil/sequences.RDS")




significantDG = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21DG.RDS")















#####FUNGIIIIIIII#######
####High abundance  High occupancy idea######
geodata  = readRDS("Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS")
ASVtab_dust<-readRDS("Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
ASVtab_dust=prune_samples(sample_data(ASVtab_dust)$point!=26, ASVtab_dust)
ASVtab_dust <- prune_taxa(taxa_sums(ASVtab_dust) > 0, ASVtab_dust)
sam_data_dust = data.frame(ASVtab_dust@sam_data) #extract dust data from phyloseq object
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "Dust_Gen"] <- "WT"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "500"] <- "Coarse SS"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "75"] <- "Medium SS"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "25"] <- "Fine SS"
sam_data_dust$Dust_Type<- factor(sam_data_dust$Dust_Type,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))




asvtab_dust = data.frame(ASVtab_dust@otu_table) #extract ASVt table not transformed






####set a threshold for the first analysis
threshold = 1






###select the ASVs based on high occupancy
selection = round((nrow(geodata)/100)*threshold) ###define the number of ASVs to pick
high_occupancy_ASVs = geodata[order(-geodata$OC_points_percent),] ### order the table
high_occupancy_ASVs = high_occupancy_ASVs[1:selection,] ### extract those with the highest occupancy
min = min(high_occupancy_ASVs$OC_points_percent) ####correct for those having the same occupancy
high_occupancy_ASVs = geodata[which(geodata$OC_points_percent >= min),] ##finally select them




Log2foldchanges = readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/log2foldchanges_fun_WT.RDS")
Log2foldchanges_occupancy  = Log2foldchanges[rownames(Log2foldchanges) %in% rownames(high_occupancy_ASVs), ]
Log2foldchanges_occupancy$Gen_spec = "Generalists"
Log2foldchanges_occupancy$abundance = high_occupancy_ASVs[match(rownames(high_occupancy_ASVs),rownames(Log2foldchanges_occupancy)),6 ]
Log2foldchanges_occupancy = head(Log2foldchanges_occupancy[order(Log2foldchanges_occupancy$abundance, decreasing=TRUE), ], 20)














###select the ASVs based on high abundance
`%nin%` = Negate(`%in%`)
geodata_specialists = geodata[rownames(geodata) %nin% rownames(high_occupancy_ASVs), ]
selection = round((nrow(geodata_specialists)/100)*threshold) ###define the number of ASVs to pick
high_abundance_ASVs = geodata_specialists[order(-geodata_specialists$Ab_bypoint),] ### order the table
high_abundance_ASVs = high_abundance_ASVs[1:selection,] ### extract those with the highest occupancy
min_ab = min(high_abundance_ASVs$Ab_bypoint) ####correct for those having the same occupancy
high_abundance_ASVs = geodata_specialists[which(geodata_specialists$Ab_bypoint >= min_ab),] ##finally select them













Log2foldchanges_abundance  = Log2foldchanges[rownames(Log2foldchanges) %in% rownames(high_abundance_ASVs), ] #select those found in high abundance
#now we need to assess how much are they present in dust samples
dust_samples_tot = colSums(asvtab_dust)
Log2foldchanges_abundance$abundance_dust = dust_samples_tot[match(rownames(Log2foldchanges_abundance),names(dust_samples_tot))] ####associate it to the ASVs
#Log2foldchanges_abundance = head(Log2foldchanges_abundance[order(Log2foldchanges_abundance$abundance_dust, decreasing=TRUE), ], 10) #select the 10 most abundant
Log2foldchanges_abundance$Gen_spec = "Specialists"
Log2foldchanges_abundance$abundance = Log2foldchanges_abundance[match(rownames(Log2foldchanges_abundance),rownames(high_abundance_ASVs)),6 ]
Log2foldchanges_abundance = head(Log2foldchanges_abundance[order(Log2foldchanges_abundance$abundance, decreasing=TRUE), ], 20)
















Log2foldchanges_abundance = Log2foldchanges_abundance[,-13]
final_log2foldchange = rbind(Log2foldchanges_occupancy, Log2foldchanges_abundance)
final_log2foldchange$dust_soil = ifelse(final_log2foldchange$log2FoldChange > 0, "Dust" , "Soil")
final_log2foldchange$name = paste(rownames(final_log2foldchange), final_log2foldchange$Genus, sep=" - ")
species = ifelse(is.na(final_log2foldchange$Species) , "sp." , final_log2foldchange$Species)
final_log2foldchange$name = paste(final_log2foldchange$name, species, sep=" ")



















fine_fungi = ggplot(final_log2foldchange, 
       aes(log2FoldChange, y = name, fill = dust_soil)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV")+ggtitle("Fungi") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_flip()+
  scale_fill_manual(name = "Sample", labels = c("Dust", "Soil"), values = c("4DBBD5B2","#996600")) + ylab(NULL)+
  facet_wrap(~ Gen_spec, scales = "free_x")
fine_fungi


ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/3rd_draft/plots/new_figures/fun_gen_spec.pdf", height = 4, width = 8)



a = ggarrange(fine_bact, fine_fungi,  common.legend = TRUE, rows = 2)
a
dev.off()



















#Associate genera to ASVs
taxa_dust = data.frame(ASVtab_dust@tax_table) #extract dust data from phyloseq object
high_occupancy_ASVs$genus = taxa_dust[match(rownames(high_occupancy_ASVs), rownames(taxa_dust)), ]
high_occupancy_ASVs$genus[is.na(high_occupancy_ASVs$genus)] <- 'unknown'
#count Nuber of unknown 
a = table(high_occupancy_ASVs$genus$Genus)
a[names(a) == "unknown"] / nrow(high_occupancy_ASVs)

high_abundance_ASVs$genus = taxa_dust[match(rownames(high_abundance_ASVs), rownames(taxa_dust)), ]
high_abundance_ASVs$genus[is.na(high_abundance_ASVs$genus)] <- 'unknown'
aggregated_ab = aggregate(Ab_total ~ genus$Genus, high_abundance_ASVs, FUN = sum)
b = table(high_abundance_ASVs$genus$Genus)
b[names(b) == "unknown"] / nrow(high_abundance_ASVs)
sequences  = readRDS("/Users/gabri/OneDrive - University of Arizona/dust_project/data/soil/sequences_fungi.RDS")
sequences[,48]

