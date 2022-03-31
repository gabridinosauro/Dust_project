library(ggplot2)
library(ggpubr)


#read sigtab
significantDG = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21DG.RDS")
#### 
increasing = significantDG[significantDG$log2FoldChange > 0, ]
###list 3 more common phyla
sort(table(increasing$Phylum) / nrow(increasing))
# Actinobacteria 35 %
# Proteobacteria 28%
# Bacteroidetes 16%
decreasing = significantDG[significantDG$log2FoldChange < 0, ]
sort(table(decreasing$Phylum) / nrow(decreasing))
# Actinobacteria 43%
# Proteobacteria 11%
# Bacteroidetes 5%
# Chloroflexi 15%

###select only the 50 most abundant ones, otherwise we cannot see anything
sigtab50 = significantDG
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name = sigtab50$Genus
sigtab50$name1 = ifelse(is.na(sigtab50$name),sigtab50$Family,sigtab50$name  )
sigtab50$name2 = ifelse(is.na(sigtab50$name1),sigtab50$Order,sigtab50$name1  )
sigtab50$name3 = ifelse(is.na(sigtab50$name2),sigtab50$Class,sigtab50$name2  )

sigtab50$name = paste(rownames(sigtab50), sigtab50$name3, sep = " - ")

# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_data.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:10,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:10,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)

DG = ggplot(increasing_decreasing, 
            aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV")+ggtitle("Bact/Arch,  WT") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_flip()+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) + ylab(NULL)
DG







significant500 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21500.RDS")
########significantDG = readRDS(file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/sigtabDEseq21DG.RDS")
###select only the 50 most abundant ones, otherwise we cannot see anything
sigtab50 = significant500
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name = sigtab50$Genus
sigtab50$name1 = ifelse(is.na(sigtab50$name),sigtab50$Family,sigtab50$name  )
sigtab50$name2 = ifelse(is.na(sigtab50$name1),sigtab50$Order,sigtab50$name1  )
sigtab50$name3 = ifelse(is.na(sigtab50$name2),sigtab50$Class,sigtab50$name2  )

sigtab50$name = paste(rownames(sigtab50), sigtab50$name3, sep = " - ")

# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_data.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:10,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:10,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)

Coarse = ggplot(increasing_decreasing, 
            aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV")+ggtitle("Bact/Arch,  coarse LP") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_flip()+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) + ylab(NULL)
Coarse


######75####
significant75 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2175.RDS")
########significantDG = readRDS(file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/sigtabDEseq21DG.RDS")
###select only the 50 most abundant ones, otherwise we cannot see anything
sigtab50 = significant75
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name = sigtab50$Genus
sigtab50$name1 = ifelse(is.na(sigtab50$name),sigtab50$Family,sigtab50$name  )
sigtab50$name2 = ifelse(is.na(sigtab50$name1),sigtab50$Order,sigtab50$name1  )
sigtab50$name3 = ifelse(is.na(sigtab50$name2),sigtab50$Class,sigtab50$name2  )

sigtab50$name = paste(rownames(sigtab50), sigtab50$name3, sep = " - ")

# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_data.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:10,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:10,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)

medium = ggplot(increasing_decreasing, 
                aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV")+ggtitle("Bact/Arch,  medium LP") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_flip()+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) + ylab(NULL)
medium



######25####
significant25 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2125.RDS")
########significantDG = readRDS(file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/sigtabDEseq21DG.RDS")
###select only the 50 most abundant ones, otherwise we cannot see anything
sigtab50 = significant25
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name = sigtab50$Genus
sigtab50$name1 = ifelse(is.na(sigtab50$name),sigtab50$Family,sigtab50$name  )
sigtab50$name2 = ifelse(is.na(sigtab50$name1),sigtab50$Order,sigtab50$name1  )
sigtab50$name3 = ifelse(is.na(sigtab50$name2),sigtab50$Class,sigtab50$name2  )

sigtab50$name = paste(rownames(sigtab50), sigtab50$name3, sep = " - ")

# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_data.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:10,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:10,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)

fine = ggplot(increasing_decreasing, 
                aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV")+ggtitle("Bact/Arch,  fine LP") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_flip()+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) + ylab(NULL)
fine

bac = ggarrange( DG, Coarse, medium, fine, common.legend = TRUE, ncol = 4)










########################################FFFFFFFFFFFFFUUUUUUUUNNNNGIIIIII##################












#read sigtab
significantDG = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21DGITS.RDS")
###select only the 50 most abundant ones, otherwise we cannot see anything
sigtab50 = significantDG
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name =  sigtab50$Genus
sigtab50$name1 = ifelse(is.na(sigtab50$name),sigtab50$Family,sigtab50$name  )
sigtab50$name2 = ifelse(is.na(sigtab50$name1),sigtab50$Order,sigtab50$name1  )
sigtab50$name3 = ifelse(is.na(sigtab50$name2),sigtab50$Class,sigtab50$name2  )
sigtab50$name4 = ifelse( sigtab50$name3 == sigtab50$Genus, paste(sigtab50$Genus, sigtab50$Species, sep = " "), paste(sigtab50$name3 ) )
sigtab50$name5 = ifelse(is.na(sigtab50$name4),sigtab50$name3,sigtab50$name4  )
sigtab50$name = paste(rownames(sigtab50), sigtab50$name4, sep = " - ")

# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:10,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:10,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)

DGfung = ggplot(increasing_decreasing, 
            aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV")+ggtitle("Fungi,  WT") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_flip()+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) + ylab(NULL)
DGfung







significant500 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21500ITS.RDS")
########significantDG = readRDS(file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/sigtabDEseq21DG.RDS")
###select only the 50 most abundant ones, otherwise we cannot see anything
sigtab50 = significant500
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name =  sigtab50$Genus
sigtab50$name1 = ifelse(is.na(sigtab50$name),sigtab50$Family,sigtab50$name  )
sigtab50$name2 = ifelse(is.na(sigtab50$name1),sigtab50$Order,sigtab50$name1  )
sigtab50$name3 = ifelse(is.na(sigtab50$name2),sigtab50$Class,sigtab50$name2  )
sigtab50$name4 = ifelse( sigtab50$name3 == sigtab50$Genus, paste(sigtab50$Genus, sigtab50$Species, sep = " "), paste(sigtab50$name3 ) )
sigtab50$name5 = ifelse(is.na(sigtab50$name4),sigtab50$name3,sigtab50$name4  )
sigtab50$name = paste(rownames(sigtab50), sigtab50$name4, sep = " - ")

# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:10,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:10,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)
increasing_decreasing = increasing_decreasing[rowSums(is.na(increasing_decreasing)) != ncol(increasing_decreasing), ]
Coarse_fung = ggplot(increasing_decreasing, 
                aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV")+ggtitle("Fungi,  coarse LP") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_flip()+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) + ylab(NULL)
Coarse_fung


######75####
significant75 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2175_ITS.RDS")
########significantDG = readRDS(file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/sigtabDEseq21DG.RDS")
###select only the 50 most abundant ones, otherwise we cannot see anything
sigtab50 = significant75
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name =  sigtab50$Genus
sigtab50$name1 = ifelse(is.na(sigtab50$name),sigtab50$Family,sigtab50$name  )
sigtab50$name2 = ifelse(is.na(sigtab50$name1),sigtab50$Order,sigtab50$name1  )
sigtab50$name3 = ifelse(is.na(sigtab50$name2),sigtab50$Class,sigtab50$name2  )
sigtab50$name4 = ifelse( sigtab50$name3 == sigtab50$Genus, paste(sigtab50$Genus, sigtab50$Species, sep = " "), paste(sigtab50$name3 ) )
sigtab50$name5 = ifelse(is.na(sigtab50$name4),sigtab50$name3,sigtab50$name4  )
sigtab50$name = paste(rownames(sigtab50), sigtab50$name4, sep = " - ")

# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:10,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:10,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)
increasing_decreasing = increasing_decreasing[rowSums(is.na(increasing_decreasing)) != ncol(increasing_decreasing), ]
medium_fung = ggplot(increasing_decreasing, 
                aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV")+ggtitle("Bact/Arch,  medium LP") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_flip()+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) + ylab(NULL)
medium_fung



######25####
significant25 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2125_ITS.RDS")
########significantDG = readRDS(file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/sigtabDEseq21DG.RDS")
###select only the 50 most abundant ones, otherwise we cannot see anything
sigtab50 = significant25
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name =  sigtab50$Genus
sigtab50$name1 = ifelse(is.na(sigtab50$name),sigtab50$Family,sigtab50$name  )
sigtab50$name2 = ifelse(is.na(sigtab50$name1),sigtab50$Order,sigtab50$name1  )
sigtab50$name3 = ifelse(is.na(sigtab50$name2),sigtab50$Class,sigtab50$name2  )
sigtab50$name4 = ifelse( sigtab50$name3 == sigtab50$Genus, paste(sigtab50$Genus, sigtab50$Species, sep = " "), paste(sigtab50$name3 ) )
sigtab50$name5 = ifelse(is.na(sigtab50$name4),sigtab50$name3,sigtab50$name4  )
sigtab50$name = paste(rownames(sigtab50), sigtab50$name4, sep = " - ")

# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:10,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:10,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)
increasing_decreasing = increasing_decreasing[rowSums(is.na(increasing_decreasing)) != ncol(increasing_decreasing), ]

fine_fung = ggplot(increasing_decreasing, 
              aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV")+ggtitle("Fungi,  fine LP") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ coord_flip()+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) + ylab(NULL)
fine_fung




fun = ggarrange(DGfung, Coarse_fung, medium_fung, fine_fung, common.legend = TRUE, ncol = 4)

#ggarrange(DG, DGfung, ncol = 1 )
# ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/diff_ab.pdf", height = 8, width = 5)


ggarrange(bac, fun, ncol = 1 )
 ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/diff_ab.pdf", height = 8, width = 10)
