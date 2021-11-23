######DEseq2 FUNGI
library(ggplot2)
library(vegan)
library(phyloseq)
library(pheatmap)
library(viridis)
library(DESeq2)
library(pheatmap)
library(metagMisc)

#########LOAD AND PREPARE THE DATA##################
asvtab<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
asvtab_no26 <- prune_taxa(taxa_sums(asvtab) > 0, asvtab)#elimiante_zero_counts
####add a column for sample repetitions (merge sampel and repetition)
sam_data = data.frame(asvtab_no26@sam_data)
sam_data = sam_data[,-1] #eliminate first column (useless)
sam_data$sample <- paste(sam_data$Point_type,sam_data$Repetition, sep = "_")
#####insert new asvtab into the phyloseq object
asvtab_no26@sam_data =sample_data(sam_data) 
asvtab_no26@sam_data$Repetition=factor(asvtab_no26@sam_data$Repetition)
rm(asvtab)
asvtab1 = data.frame(t(asvtab_no26@otu_table))
asvtab1 = asvtab1[rowSums(asvtab1 != 0) >= 3, ]
asvtab_no26@otu_table =otu_table(asvtab1, taxa_are_rows = TRUE)




#####RUN DESeq2##############################################################
diagdds = phyloseq_to_deseq2(asvtab_no26, ~  sample +  Dust_Type) #convert phyloseq into deseq2 ##samples are paired!
diagdds$Dust_Type <- relevel(diagdds$Dust_Type, ref = "Soil") ####relevel placing soil as reference

#filtering (remove certain sequences)
diagdds = estimateSizeFactors(diagdds) #prestep for filtering
##### equal or more than 3 samples (last number)samples with normalized counts higher than 2 (first number)
idx <- rowSums( counts(diagdds, normalized=TRUE) >= 2 ) >=3  
diagdds <- diagdds[idx,] 
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
resultsNames(diagdds)
#saveRDS(diagdds, file = "/Users/gabri/OneDrive - University of Arizona/dust_project/data/dust/DEseq2file.RDS")


#############################RESULTS DG (dust_generator) fungi#############################
res = results(diagdds,  contrast=c("Dust_Type","Dust_Gen","Soil"), pAdjustMethod = "BH", alpha =  0.05)
summary(res)
resLFC <- lfcShrink(diagdds, coef=c("Dust_Type_Dust_Gen_vs_Soil"), type="apeglm")
summary(resLFC, pAdjustMethod = "BH", alpha =  0.05) #85 up , 103 down
par(mfrow=c(1,1))
DESeq2::plotMA(resLFC)
Log2foldchanges = data.frame(resLFC)
Log2foldchanges = cbind(Log2foldchanges, as(tax_table(asvtab_no26)[rownames(Log2foldchanges), ], "matrix"))
saveRDS(Log2foldchanges, "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/log2foldchanges_fun_WT.RDS")


ggplot(as(resLFC, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtabDG = cbind(as(sigtab, "data.frame"), as(tax_table(asvtab_no26)[rownames(sigtab), ], "matrix"))
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
for ( i in 10:18) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ###tarting problems ( e.g. 1 read in soil...)
for ( i in 19:27) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ###problems ( e.g. 1 read in soil...)
for ( i in 26:34) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ######problems continuing (e.g. 1 read in soil)

#####Let us see what happens in the increasing ones
a = a1[a1$log2FoldChange > 0,]
a = a[order(-a$baseMean),]
par(mfrow=c(3,3))
for ( i in 1:9) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ##till 9, solid
for ( i in 10:18) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ### till 18 ok
for ( i in 19:27) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ###ok
for ( i in 26:34) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ######ok
for ( i in 35:43) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ###starting to show things with one count
for ( i in 44:52) {plotCounts(diagdds, gene=rownames(a)[i], intgroup="Dust_Type") } ###problems (1 count)

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
sum(sigtabDG$log2FoldChange>0) #41
sum(sigtabDG$log2FoldChange<0) #14
saveRDS(sigtabDG, file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21DGITS.RDS")


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
col1 = metadata$Point_type
col1[col1=="Santa_Cruz_River_bank"]<-"SCR"
col1[col1=="Picacho_Peak"]<-"PP"
col1[col1=="Desert_shrub_2"]<-"DS2"
col1[col1=="Desert_shrub_1"]<-"DS1"

# Data frame with column annotations.
mat_col <- data.frame(Sample = col)

#rownames(mat_col) <- colnames(asvtab5)
mat_colors <- list(Sample = c("blue","red"), Sampling_point = c( PP = "#D16103" ,SCR ="#00AFBB" , DS2 = "#E7B800",DS1 ="#52854C"))
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
  main              = "Fungi, WT samples",
  cluster_cols =  F,
  treeheight_row = 0
)
a



#############################RESULTS 500 fungi#############################

res = results(diagdds,  contrast=c("Dust_Type","500","Soil"), pAdjustMethod = "BH", alpha =  0.05)
summary(res)
resLFC <- lfcShrink(diagdds, coef=c("Dust_Type_500_vs_Soil"), type="apeglm")
summary(resLFC, pAdjustMethod = "BH", alpha =  0.05) #85 up , 103 down
par(mfrow=c(1,1))
DESeq2::plotMA(resLFC)
ggplot(as(resLFC, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtab500 = cbind(as(sigtab, "data.frame"), as(tax_table(asvtab_no26)[rownames(sigtab), ], "matrix"))
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
sigtab500 = data.frame(sigtab500[rownames(sigtab500) %in% soil_guys |rownames(sigtab500) %in% dust_guys, ])
sum(sigtab500$log2FoldChange>0) #3
sum(sigtab500$log2FoldChange<0) #9
saveRDS(sigtab500, file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21500ITS.RDS")

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
mat_colors <- list(Sample = c("blue","red"), Sampling_point = c( PP = "#D16103" ,SCR ="#00AFBB" , DS2 = "#E7B800",DS1 ="#52854C"))
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
Species = sigtab500[match(ASV_names, rownames(sigtab500)),12]
Genus = sigtab500[match(ASV_names, rownames(sigtab500)),11]
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
  #fontsize          = 14,
  main              = "Fungi, Coarse SS samples",
  cluster_cols =  F,
  treeheight_row = 0
)
a






#############################RESULTS 75 fungi#############################

res = results(diagdds,  contrast=c("Dust_Type","75","Soil"), pAdjustMethod = "BH", alpha =  0.05)
summary(res)
resLFC <- lfcShrink(diagdds, coef=c("Dust_Type_75_vs_Soil"), type="apeglm")
summary(resLFC, pAdjustMethod = "BH", alpha =  0.05) #85 up , 103 down
par(mfrow=c(1,1))
DESeq2::plotMA(resLFC)
ggplot(as(resLFC, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtab75 = cbind(as(sigtab, "data.frame"), as(tax_table(asvtab_no26)[rownames(sigtab), ], "matrix"))
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
sum(sigtab75$log2FoldChange>0) #5
sum(sigtab75$log2FoldChange<0) #10
saveRDS(sigtab75, file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2175_ITS.RDS")


#######################heatmap 75################################

asvtab2 = asvtab1[rownames(asvtab1) %in% rownames(sigtab75),] #####select only significant ones
asvtab3 <- sweep(asvtab2, 2, colSums(asvtab1), FUN = "/")
asvtab3 <- as.data.frame(t(asvtab3))
asvtab4 <- decostand(asvtab3, "standardize")
asvtab4 = asvtab4
asvtab5 = t(asvtab4)



#########match it with metadata%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colnames(asvtab5) = rownames(sam_data)
metadata = sam_data[which(sam_data$Dust_Type == "Soil" | sam_data$Dust_Type == "75"),]
mat = as.matrix(asvtab5[,colnames(asvtab5) %in% rownames(metadata)])

###############Change the names of the sampling points- dust%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col = metadata$Dust_Type
levels(col)[levels(col)=="75"] <- "Medium SS"
col1 = metadata$Point_type
col1[col1=="Santa_Cruz_River_bank"]<-"SCR"
col1[col1=="Picacho_Peak"]<-"PP"
col1[col1=="Desert_shrub_2"]<-"DS2"
col1[col1=="Desert_shrub_1"]<-"DS1"

# Data frame with column annotations.
mat_col <- data.frame(Sample = col)

#rownames(mat_col) <- colnames(asvtab5)
mat_colors <- list(Sample = c("blue","red"), Sampling_point = c( PP = "#D16103" ,SCR ="#00AFBB" , DS2 = "#E7B800",DS1 ="#52854C"))
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
Species = sigtab75[match(ASV_names, rownames(sigtab75)),13]
Genus = sigtab75[match(ASV_names, rownames(sigtab75)),12]
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
  #fontsize          = 14,
  main              = "Fungi, Medium SS samples",
  cluster_cols =  F,
  treeheight_row = 0
)
a






#############################RESULTS 25 fungi#############################

res = results(diagdds,  contrast=c("Dust_Type","25","Soil"), pAdjustMethod = "BH", alpha =  0.05)
summary(res)
resLFC <- lfcShrink(diagdds, coef=c("Dust_Type_25_vs_Soil"), type="apeglm")
summary(resLFC, pAdjustMethod = "BH", alpha =  0.05) #85 up , 103 down
par(mfrow=c(1,1))
DESeq2::plotMA(resLFC)
ggplot(as(resLFC, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtab25 = cbind(as(sigtab, "data.frame"), as(tax_table(asvtab_no26)[rownames(sigtab), ], "matrix"))
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
sum(sigtab25$log2FoldChange>0) #9
sum(sigtab25$log2FoldChange<0) #11
saveRDS(sigtab25, file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2125_ITS.RDS")


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
mat_colors <- list(Sample = c("blue","red"), Sampling_point = c( PP = "#D16103" ,SCR ="#00AFBB" , DS2 = "#E7B800",DS1 ="#52854C"))
col[col=="25"]<-"Fine SS"
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
Species = sigtab25[match(ASV_names, rownames(sigtab25)),13]
Genus = sigtab25[match(ASV_names, rownames(sigtab25)),12]
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
  #fontsize          = 14,
  main              = "Fungi, Fine SS samples",
  cluster_cols =  F,
  treeheight_row = 0
)
a


