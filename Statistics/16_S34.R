####DEseq2 Genus level####
#######DEseq2 Analysis dust project Gabri ####diagnostic script
library(ggplot2)
library(vegan)
library(phyloseq)
library(pheatmap)
library(viridis)
library(DESeq2)
library(pheatmap)
library(metagMisc)








#########LOAD AND PREPARE THE DATA##################
ASVtab_dust<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
ASVtab_dust_gen = tax_glom(ASVtab_dust, taxrank=rank_names(ASVtab_dust)[6], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
#####insert new asvtab into the phyloseq object
sam_data = data.frame(ASVtab_dust_gen@sam_data) #extract sam_data
sam_data = sam_data[,-1] #eliminate first column (useless)
sam_data$sample <- paste(sam_data$Point_type,sam_data$Repetition, sep = "_")##create a repetition column
#####insert new asvtab into the phyloseq object
ASVtab_dust_gen@sam_data =sample_data(sam_data) 
ASVtab_dust_gen@sam_data$Repetition=factor(ASVtab_dust_gen@sam_data$Repetition)
asvtab1 = data.frame(t(ASVtab_dust_gen@otu_table)) # extract ASV table for later












#####RUN DESeq2##############################################################
diagdds = phyloseq_to_deseq2(ASVtab_dust_gen, ~  sample +  Dust_Type) #convert phyloseq into deseq2 ##samples are paired!
diagdds$Dust_Type <- relevel(diagdds$Dust_Type, ref = "Soil") ####relevel placing soil as reference
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
resultsNames(diagdds)






#############################RESULTS DG (dust_generator) bacteria#############################
res = results(diagdds,  name=c("Dust_Type_Dust_Gen_vs_Soil"), pAdjustMethod = "BH", alpha =  0.05)
summary(res)
taxa_tab = data.frame(ASVtab_dust@tax_table)
resLFC <- lfcShrink(diagdds, coef=c("Dust_Type_Dust_Gen_vs_Soil"), type="apeglm")
res_gens = data.frame(resLFC)
rownames(res_gens) = taxa_tab$Genus[match(rownames(res_gens), rownames(taxa_tab))]
res_gens = res_gens[res_gens$padj<0.05,]

res_gens = res_gens[1:50,]
#write.table(res_gens, "/Users/gabri/OneDrive - University of Arizona/dust_project/table_genus.txt")


resLFC <- lfcShrink(diagdds, coef=c("Dust_Type_Dust_Gen_vs_Soil"), type="apeglm")
summary(resLFC, pAdjustMethod = "BH", alpha =  0.1) #292 up , 289 down
par(mfrow=c(1,1))
DESeq2::plotMA(resLFC)
ggplot(as(resLFC, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.05, fill = "Royalblue", boundary = 0)
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtabDG = cbind(as(sigtab, "data.frame"), as(tax_table(ASVtab_dust_gen)[rownames(sigtab), ], "matrix"))
head(sigtabDG)
dim(sigtabDG)




###############DIAGNOSTIC PLOTS##########################
######## 9 less abundant
par(mfrow=c(3,3))
a1 = sigtabDG[order(sigtabDG$baseMean),]
for ( i in 1:9) {plotCounts(diagdds, gene=rownames(a1)[i], intgroup="Dust_Type") }
####plot 9 most abundant ASVs
par(mfrow=c(3,3))
a1 = sigtabDG[order(-sigtabDG$baseMean),]
for ( i in 1:9) {plotCounts(diagdds, gene=rownames(a1)[i], intgroup="Dust_Type") }
####Plot only the negative ones
a = a1[a1$log2FoldChange < 0,]
a = a[order(-a$baseMean),]
par(mfrow=c(3,3))
for ( i in 1:9) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ##till 9, solid
for ( i in 10:18) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ### solid
for ( i in 19:27) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ###solid
for ( i in 26:34) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ######starting problems
for ( i in 35:43) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ########## problems continuing (e.g. 1 read in soil)
for ( i in 44:52) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ########problems, bad data

#####Let us see what happens in the increasing ones
a = a1[a1$log2FoldChange > 0,]
a = a[order(-a$baseMean),]
par(mfrow=c(3,3))
for ( i in 1:9) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ##till 9, solid
for ( i in 10:18) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ### till 18 ok
for ( i in 19:27) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ###ok
for ( i in 26:34) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ######ok
for ( i in 35:43) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ###prblems
for ( i in 44:52) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ###problems

######################################Filtering############################################################
####I will filter the results, only considering those ASVs, which have more than 2 samples in soil (decreasing, and in dust (increasing))

###soil filtering
database = counts(diagdds, normalized = TRUE)#extract dataframe from DEseq object
database = database[, colnames(database) %in% c(rownames(sam_data[sam_data$Dust_Type == "Soil",]))]
soil_guys = rownames(database[rowSums(database != 0) >= 3, ])
####dust_filtering
database = counts(diagdds, normalized = TRUE)#extract dataframe from DEseq object
database = database[, colnames(database) %in% c(rownames(sam_data[sam_data$Dust_Type == "Dust_Gen",]))]
dust_guys = rownames(database[rowSums(database != 0) >= 3, ])
rm(database)


###save new filter out from sigtab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigtabDG = sigtabDG[rownames(sigtabDG) %in% soil_guys |rownames(sigtabDG) %in% dust_guys, ]
sum(sigtabDG$log2FoldChange>0) #48
sum(sigtabDG$log2FoldChange<0) # 32
#saveRDS(sigtabDG, file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/DEseq2_bac_DG_genus.RDS")


############HEATMAP#######################
############Relative abundances and normalization
asvtab2 = asvtab1[rownames(asvtab1) %in% rownames(sigtabDG),] #####select only significant ones
asvtab3 <- sweep(asvtab2, 2, colSums(asvtab1), FUN = "/")
asvtab3 <- as.data.frame(t(asvtab3))
asvtab4 <- decostand(asvtab3, "standardize")
asvtab4 = asvtab4
asvtab5 = t(asvtab4)



#########match it with metadata%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colnames(asvtab5) = rownames(sam_data)
metadata = sam_data[which(sam_data$Dust_Type == "Soil" | sam_data$Dust_Type == "Dust_Gen"),]
mat = as.matrix(asvtab5[,colnames(asvtab5) %in% rownames(metadata)])

###############Change the names of the sampling points- dust%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col = metadata$Dust_Type
levels(col)[levels(col)=="Dust_Gen"] <- "WT"
#col[col=="Dust_Gen"]<-as.factor("Dust")
col1 = metadata$Point_type
col1[col1=="Santa_Cruz_River_bank"]<-"SCR"
col1[col1=="Picacho_Peak"]<-"PP"
col1[col1=="Desert_shrub_2"]<-"DS2"
col1[col1=="Desert_shrub_1"]<-"DS1"

# Data frame with column annotations.
mat_col <- data.frame(Sample = col)

#rownames(mat_col) <- colnames(asvtab5)
mat_colors <- list(Sample = c("blue","red"), Sampling_point = c( PP = "#D16103" ,SCR ="#00AFBB" , DS1 = "#E7B800",DS2 ="#52854C"))
col[col=="Dust_Gen"]<-"WT"
names(mat_colors$Sample) <- unique(col)

######plot it%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#mat <- mat[, order(col[1], decreasing = TRUE)] 
# Data frame with column annotations.
mat_col <- data.frame(Sample = col, Sampling_point = col1)
rownames(mat_col) <- colnames(mat)


soil_dust = metadata[order(metadata$Dust_Type),]
soil_dust = rownames(soil_dust)
newmat <- mat[, match(soil_dust, colnames(mat)) ]        # Reorder data frame

###Change the rownames into species /genus
ASV_names = rownames(newmat)
Species = sigtabDG[match(ASV_names, rownames(sigtabDG)),12]
Genus = sigtabDG[match(ASV_names, rownames(sigtabDG)),11]
full_names = paste(ASV_names, Genus, Species, sep =' ')
rownames(newmat) = full_names


range <- max(abs(newmat))
#pheatmap(newmat, breaks = seq(-range/3, range, length.out = 100))
#tiff("/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/test3.tiff", units="in", width=7, height=10, res=300)

a = pheatmap(
  mat               = newmat,
  #color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 4,
  main              = "Differentially Expressed ASVs, Bacteria and Archaea, WT samples",
  cluster_cols =  F,
  treeheight_row = 0
)

#a

#dev.off()
######plot the increase decrease patterns of the most abundant ones
###order the table by abundance

####select only the 50 most abundant ones, otherwise we cannot see anything
sigtab50 = sigtabDG
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name = paste( sigtab50$Genus, "spp", sep=" ")

###plot it
ggplot(sigtab50, 
       aes(log2FoldChange, y = name, fill = color)) + 
  geom_col(border = NA) + theme_bw() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV") +ggtitle("Bact/Arch, differentially abundant Genera, WT") +
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2"))
dev.off()
# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/geo_dist_data.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:5,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:5,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)

DG = ggplot(increasing_decreasing, 
       aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("Genus")+ggtitle("Bact/Arch,  WT") + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
   scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) 
DG




#############################RESULTS 500 bacteria#############################

res = results(diagdds,  contrast=c("Dust_Type","500","Soil"), pAdjustMethod = "BH", alpha =  0.1)
summary(res) #103 up , 106 down
par(mfrow=c(1,1))
resLFC <- lfcShrink(diagdds, coef=c("Dust_Type_500_vs_Soil"), type="apeglm")
DESeq2::plotMA(resLFC)
ggplot(as(resLFC, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
alpha = 0.05
sigtab = res[which(resLFC$padj < alpha), ]
sigtab500 = cbind(as(sigtab, "data.frame"), as(tax_table(ASVtab_dust_gen)[rownames(sigtab), ], "matrix"))
head(sigtab500)
dim(sigtab500)




###soil filtering
database = counts(diagdds, normalized = TRUE)#extract dataframe from DEseq object
database = database[, colnames(database) %in% c(rownames(sam_data[sam_data$Dust_Type == "Soil",]))]
soil_guys = rownames(database[rowSums(database != 0) >= 3, ])
####dust_filtering
database = counts(diagdds, normalized = TRUE)#extract dataframe from DEseq object
database = database[, colnames(database) %in% c(rownames(sam_data[sam_data$Dust_Type == "500",]))]
dust_guys = rownames(database[rowSums(database != 0) >= 3, ])
rm(database)

###save new filter out from sigtab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigtab500 = sigtab500[rownames(sigtab500) %in% soil_guys |rownames(sigtab500) %in% dust_guys, ]
sum(sigtab500$log2FoldChange>0) #23
sum(sigtab500$log2FoldChange<0) #17
#saveRDS(sigtab500, file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/sigtabDEseq21500.RDS")

##########heatmap 500#############

asvtab2 = asvtab1[rownames(asvtab1) %in% rownames(sigtab500),] #####select only significant ones
asvtab3 <- sweep(asvtab2, 2, colSums(asvtab1), FUN = "/")
asvtab3 <- as.data.frame(t(asvtab3))
asvtab4 <- decostand(asvtab3, "standardize")
asvtab4 = asvtab4
asvtab5 = t(asvtab4)



#########match it with metadata%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colnames(asvtab5) = rownames(sam_data)
metadata = sam_data[which(sam_data$Dust_Type == "Soil" | sam_data$Dust_Type == "500"),]
mat = as.matrix(asvtab5[,colnames(asvtab5) %in% rownames(metadata)])

###############Change the names of the sampling points- dust%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col = metadata$Dust_Type
levels(col)[levels(col)=="500"] <- "Coarse SS"
col1 = metadata$Point_type
col1[col1=="Santa_Cruz_River_bank"]<-"SCR"
col1[col1=="Picacho_Peak"]<-"PP"
col1[col1=="Desert_shrub_2"]<-"DS2"
col1[col1=="Desert_shrub_1"]<-"DS1"

# Data frame with column annotations.
mat_col <- data.frame(Sample = col)

#rownames(mat_col) <- colnames(asvtab5)
mat_colors <- list(Sample = c("blue","red"), Sampling_point = c( PP = "#D16103" ,SCR ="#00AFBB" , DS1 = "#E7B800",DS2 ="#52854C"))
col[col=="500"]<-"Coarse SS"
names(mat_colors$Sample) <- unique(col)

######plot it%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#mat <- mat[, order(col[1], decreasing = TRUE)] 
# Data frame with column annotations.
mat_col <- data.frame(Sample = col, Sampling_point = col1)
rownames(mat_col) <- colnames(mat)


soil_dust = metadata[order(metadata$Dust_Type),]
soil_dust = rownames(soil_dust)
newmat <- mat[, match(soil_dust, colnames(mat)) ]        # Reorder data frame


###Change the rownames into species /genus
ASV_names = rownames(newmat)
Species = sigtab500[match(ASV_names, rownames(sigtab500)),13]
Genus = sigtab500[match(ASV_names, rownames(sigtab500)),12]
full_names = paste(ASV_names, Genus, Species, sep =' ')
rownames(newmat) = full_names




range <- max(abs(newmat))
#pheatmap(newmat, breaks = seq(-range/3, range, length.out = 100))
#tiff("/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/test6.tiff", units="in", width=8.5, height=10, res=300)


a = pheatmap(
  mat               = newmat,
  #color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  #fontsize          = 14,
  main              = "Bacteria and Archaea, Coarse SS samples",
  cluster_cols =  F,
  treeheight_row = 0
)
#a


dev.off()

####select only the 50 most abundant ones, otherwise we cannot see anything
sigtab50 = sigtab500
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name = paste( sigtab50$Genus, "spp", sep=" ")

###plot it
ggplot(sigtab50, 
       aes(log2FoldChange, y = name, fill = color)) + 
  geom_col(border = NA) + theme_bw() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV") +ggtitle("Bact/Arch, differentially abundant Genera, Coarse SS")
dev.off()
# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/geo_dist_data.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:5,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:5,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)

coarse = ggplot(increasing_decreasing, 
            aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("Genus")+ggtitle("Bact/Arch,  coarse SS") + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) 
coarse



#############################RESULTS 75 bacteria#############################

res = results(diagdds,  contrast=c("Dust_Type","75","Soil"), pAdjustMethod = "BH", alpha =  0.1)
summary(res) #114 up , 102 down
par(mfrow=c(1,1))
resLFC <- lfcShrink(diagdds, coef=c("Dust_Type_75_vs_Soil"), type="apeglm")
DESeq2::plotMA(resLFC)
ggplot(as(res, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtab75 = cbind(as(sigtab, "data.frame"), as(tax_table(ASVtab_dust_gen)[rownames(sigtab), ], "matrix"))
head(sigtab75)
dim(sigtab75)




###soil filtering
database = counts(diagdds, normalized = TRUE)#extract dataframe from DEseq object
database = database[, colnames(database) %in% c(rownames(sam_data[sam_data$Dust_Type == "Soil",]))]
soil_guys = rownames(database[rowSums(database != 0) >= 3, ])
####dust_filtering
database = counts(diagdds, normalized = TRUE)#extract dataframe from DEseq object
database = database[, colnames(database) %in% c(rownames(sam_data[sam_data$Dust_Type == "75",]))]
dust_guys = rownames(database[rowSums(database != 0) >= 3, ])
rm(database)

###save new filter out from sigtab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigtab75 = sigtab75[rownames(sigtab75) %in% soil_guys |rownames(sigtab75) %in% dust_guys, ]
sum(sigtab75$log2FoldChange>0) #28
sum(sigtab75$log2FoldChange<0) #33
#saveRDS(sigtab75, file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/sigtabDEseq2175.RDS")


#########match it with metadata%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colnames(asvtab5) = rownames(sam_data)
metadata = sam_data[which(sam_data$Dust_Type == "Soil" | sam_data$Dust_Type == "75"),]
mat = as.matrix(asvtab5[,colnames(asvtab5) %in% rownames(metadata)])

###############Change the names of the sampling points- dust%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col = metadata$Dust_Type
levels(col)[levels(col)=="75"] <- "Medium SS"
#col[col=="Dust_Gen"]<-as.factor("Dust")
col1 = metadata$Point_type
col1[col1=="Santa_Cruz_River_bank"]<-"SCR"
col1[col1=="Picacho_Peak"]<-"PP"
col1[col1=="Desert_shrub_2"]<-"DS2"
col1[col1=="Desert_shrub_1"]<-"DS1"

# Data frame with column annotations.
mat_col <- data.frame(Sample = col)

#rownames(mat_col) <- colnames(asvtab5)
mat_colors <- list(Sample = c("blue","red"), Sampling_point = c( PP = "#D16103" ,SCR ="#00AFBB" , DS1 = "#E7B800",DS2 ="#52854C"))
col[col=="75"]<-"Medium SS"
names(mat_colors$Sample) <- unique(col)

######plot it%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#mat <- mat[, order(col[1], decreasing = TRUE)] 
# Data frame with column annotations.
mat_col <- data.frame(Sample = col, Sampling_point = col1)
rownames(mat_col) <- colnames(mat)


soil_dust = metadata[order(metadata$Dust_Type),]
soil_dust = rownames(soil_dust)
newmat <- mat[, match(soil_dust, colnames(mat)) ]        # Reorder data frame



###Change the rownames into species /genus
ASV_names = rownames(newmat)
Species = sigtab75[match(ASV_names, rownames(sigtab75)),12]
Genus = sigtab75[match(ASV_names, rownames(sigtab75)),11]
full_names = paste(ASV_names, Genus, Species, sep =' ')
rownames(newmat) = full_names






range <- max(abs(newmat))
#pheatmap(newmat, breaks = seq(-range/3, range, length.out = 100))

#tiff("/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/test7.tiff", units="in", width=8.5, height=10, res=300)

a = pheatmap(
  mat               = newmat,
  #color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 6,
  main              = "Bacteria and Archaea, Medium SS samples",
  cluster_cols =  F,
  treeheight_row = 0
)

a
dev.off()

sigtab50 = sigtab75
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name = paste( sigtab50$Genus, "spp", sep=" ")

###plot it
ggplot(sigtab50, 
       aes(log2FoldChange, y = name, fill = color)) + 
  geom_col(border = NA) + theme_bw() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV") +ggtitle("Bact/Arch, differentially abundant Genera, Medium SP")
dev.off()
# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/geo_dist_data.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:5,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:5,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)

Medium = ggplot(increasing_decreasing, 
            aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("Genus")+ggtitle("Bact/Arch,  medium SS") + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) 
Medium

#############################RESULTS 25 bacteria#############################

res = results(diagdds,  contrast=c("Dust_Type","25","Soil"), pAdjustMethod = "BH", alpha =  0.1)
resLFC <- lfcShrink(diagdds, coef=c("Dust_Type_25_vs_Soil"), type="apeglm")
summary(resLFC) #36 # 54
par(mfrow=c(1,1))
DESeq2::plotMA(resLFC)
ggplot(as(resLFC, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtab25 = cbind(as(sigtab, "data.frame"), as(tax_table(ASVtab_dust_gen)[rownames(sigtab), ], "matrix"))
head(sigtab25)
dim(sigtab25)




###soil filtering
database = counts(diagdds, normalized = TRUE)#extract dataframe from DEseq object
database = database[, colnames(database) %in% c(rownames(sam_data[sam_data$Dust_Type == "Soil",]))]
soil_guys = rownames(database[rowSums(database != 0) >= 3, ])
####dust_filtering
database = counts(diagdds, normalized = TRUE)#extract dataframe from DEseq object
database = database[, colnames(database) %in% c(rownames(sam_data[sam_data$Dust_Type == "25",]))]
dust_guys = rownames(database[rowSums(database != 0) >= 3, ])
rm(database)

###save new filter out from sigtab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigtab25 = sigtab25[rownames(sigtab25) %in% soil_guys |rownames(sigtab25) %in% dust_guys, ]
sum(sigtab25$log2FoldChange>0) #67
sum(sigtab25$log2FoldChange<0) #46
#saveRDS(sigtab25, file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/sigtabDEseq2125.RDS")


###############heatmap 25#########
asvtab2 = asvtab1[rownames(asvtab1) %in% rownames(sigtab25),] #####select only significant ones
asvtab3 <- sweep(asvtab2, 2, colSums(asvtab1), FUN = "/")
asvtab3 <- as.data.frame(t(asvtab3))
asvtab4 <- decostand(asvtab3, "standardize")
asvtab4 = asvtab4
asvtab5 = t(asvtab4)



#########match it with metadata%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colnames(asvtab5) = rownames(sam_data)
metadata = sam_data[which(sam_data$Dust_Type == "Soil" | sam_data$Dust_Type == "25"),]
mat = as.matrix(asvtab5[,colnames(asvtab5) %in% rownames(metadata)])

###############Change the names of the sampling points- dust%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col = metadata$Dust_Type
levels(col)[levels(col)=="25"] <- "Fine SS"
col1 = metadata$Point_type
col1[col1=="Santa_Cruz_River_bank"]<-"SCR"
col1[col1=="Picacho_Peak"]<-"PP"
col1[col1=="Desert_shrub_2"]<-"DS2"
col1[col1=="Desert_shrub_1"]<-"DS1"

# Data frame with column annotations.
mat_col <- data.frame(Sample = col)

#rownames(mat_col) <- colnames(asvtab5)
mat_colors <- list(Sample = c("blue","red"), Sampling_point = c( PP = "#D16103" ,SCR ="#00AFBB" , DS1 = "#E7B800",DS2 ="#52854C"))
col[col=="25"]<-"Fine SP"
names(mat_colors$Sample) <- unique(col)

######plot it%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#mat <- mat[, order(col[1], decreasing = TRUE)] 
# Data frame with column annotations.
mat_col <- data.frame(Sample = col, Sampling_point = col1)
rownames(mat_col) <- colnames(mat)


soil_dust = metadata[order(metadata$Dust_Type),]
soil_dust = rownames(soil_dust)
newmat <- mat[, match(soil_dust, colnames(mat)) ]        # Reorder data frame


###Change the rownames into species /genus
ASV_names = rownames(newmat)
Species = sigtab25[match(ASV_names, rownames(sigtab25)),12]
Genus = sigtab25[match(ASV_names, rownames(sigtab25)),11]
full_names = paste(ASV_names, Genus, Species, sep =' ')
rownames(newmat) = full_names




range <- max(abs(newmat))
#pheatmap(newmat, breaks = seq(-range/3, range, length.out = 100))


a = pheatmap(
  mat               = newmat,
  #color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 8,
  main              = "Bacteria and Archaea, Fine SS samples",
  cluster_cols =  F
)
a
dev.off()


sigtab50 = sigtab25
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab50$name = paste( sigtab50$Genus, "spp", sep=" ")

###plot it
ggplot(sigtab50, 
       aes(log2FoldChange, y = name, fill = color)) + 
  geom_col(border = NA) + theme_bw() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("ASV") +ggtitle("Bact/Arch, differentially abundant Genera, Fine SP")
dev.off()
# + geom_errorbar( aes(xmin=log2FoldChange-lfcSE , xmax=log2FoldChange+lfcSE , y = name), width=0.4, colour="orange", alpha=0.9, size=0.5)

geodata =  readRDS("/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/geo_dist_data.RDS")
sigtab50$abundance = geodata$Ab_total[match(row.names(sigtab50), row.names(geodata))] ### match the dataframe

#####Increasing-decreasing
sigtab50$color = ifelse(sigtab50$log2FoldChange > 0, "Increase" , "Decrease")
increasing = sigtab50[sigtab50$color == "Increase", ]
decreasing = sigtab50[sigtab50$color == "Decrease", ]

###Order them by abundance and select first ten
increasing = increasing[sort(-increasing$abundance),]
increasing = increasing[1:5,]
decreasing = decreasing[sort(-decreasing$abundance),]
decreasing = decreasing[1:5,]
###bind them
increasing_decreasing = rbind(increasing,decreasing)

Fine = ggplot(increasing_decreasing, 
                aes(log2FoldChange, y = name, fill = color)) +
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3) +
  ylab("Genus")+ggtitle("Bact/Arch,  fine SS") + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","4DBBD5B2")) 
Fine


ggarrange(DG, coarse, Medium, Fine, common.legend = TRUE)





