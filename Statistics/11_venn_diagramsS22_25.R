library(ggvenn)



####Common ASVs between conditions#######
significant500 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21500.RDS")
significant75 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2175.RDS")
significant25 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2125.RDS")
significantDG = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21DG.RDS")










####Increasing#######
increasing500 = rownames(significant500[which(significant500$log2FoldChange > 0 ),]) 
increasing75  = rownames(significant75[which(significant75$log2FoldChange > 0 ),])
increasing25  = rownames(significant25[which(significant25$log2FoldChange > 0 ),])
increasingDG   = rownames(significantDG[which(significantDG$log2FoldChange > 0 ),])







####ASV found increasing in all dust conditions
Increased_shared_all = Reduce(intersect, list(increasing500,increasing75,increasing25,increasingDG ) )#none
Increased_shared_all 
#[1] "ASV13"   "ASV57"   "ASV65"   "ASV101"  "ASV138"  "ASV156"  "ASV285"  "ASV554"  "ASV1699"











#####Venn Diagram
list_increasing = list(`Coarse LP` =increasing500, `Medium LP` = increasing75, `Fine LP` =  increasing25, `WT` = increasingDG )
a = ggvenn(
  list_increasing, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
a












library(gplots)
ItemsList <- venn(list_increasing, show.plot = FALSE)
attributes(ItemsList)$intersections$WT
##Increasing in all dust samples
shared_all= significantDG[rownames(significantDG) %in% Increased_shared_all,]
###Increasing only in WT (10 most abundant)
Increased_WT = significantDG[rownames(significantDG) %in% attributes(ItemsList)$intersections$WT,]
Increased_WT =   Increased_WT[order(-Increased_WT$baseMean),] 
Increased_WT = Increased_WT[1:10,]
###Increasing only in 500 (10 most abundant)
Increased_500 = significant500[rownames(significant500) %in% attributes(ItemsList)$intersections$`Coarse SS`,]
Increased_500 =   Increased_500[order(-Increased_500$baseMean),] 
Increased_500 = Increased_500[1:10,]
###Increasing only in 75 (10 most abundant)
Increased_75 = significant75[rownames(significant75) %in% attributes(ItemsList)$intersections$`Medium SS`,]
Increased_75 =   Increased_75[order(-Increased_75$baseMean),] 
Increased_75 = Increased_75[1:10,]
###Increasing only in 25 (10 most abundant)
Increased_25 = significant25[rownames(significant25) %in% attributes(ItemsList)$intersections$`Fine SS`,]
Increased_25 =   Increased_25[order(-Increased_25$baseMean),] 
Increased_25 = Increased_25[1:10,]




















####Decreasing######
decreasing500 = rownames(significant500[which(significant500$log2FoldChange < 0 ),]) 
decreasing75  = rownames(significant75[which(significant75$log2FoldChange < 0 ),])
decreasing25  = rownames(significant25[which(significant25$log2FoldChange < 0 ),])
decreasingDG   = rownames(significantDG[which(significantDG$log2FoldChange < 0 ),])
decreased_shared_all = Reduce(intersect, list(decreasing500,decreasing75,decreasing25,decreasingDG ) )#none
decreased_shared_all # "ASV128"  "ASV163"  "ASV300"  "ASV396"  "ASV474"  "ASV853"  "ASV2374" "ASV3535" "ASV3565" "ASV3968" "ASV7649" "ASV8151"









#####Venn Diagram
list_decreasing = list(`Coarse LP` =decreasing500, `Medium LP` = decreasing75, `Fine LP` =  decreasing25, `WT` = decreasingDG )
b  = ggvenn(
  list_decreasing, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
ggarrange(a,b, labels = "AUTO")






















############FUNGI##############
###Common ASVs between conditions#######

significant500 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21500ITS.RDS")
significant75 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2175_ITS.RDS")
significant25 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2125_ITS.RDS")
significantDG = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21DGITS.RDS")








####Increasing#######
increasing500 = rownames(significant500[which(significant500$log2FoldChange > 0 ),]) 
increasing75  = rownames(significant75[which(significant75$log2FoldChange > 0 ),])
increasing25  = rownames(significant25[which(significant25$log2FoldChange > 0 ),])
increasingDG   = rownames(significantDG[which(significantDG$log2FoldChange > 0 ),])











####ASV found increasing in all dust conditions
Increased_shared_all = Reduce(intersect, list(increasing500,increasing75,increasing25,increasingDG ) )#none
Increased_shared_all 
#[1] "ASVfun266"






#####Venn Diagram
list_increasing = list(`Coarse LP` =increasing500, `Medium LP` = increasing75, `Fine LP` =  increasing25, `WT` = increasingDG )
a = ggvenn(
  list_increasing, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
a
library(gplots)
ItemsList <- venn(list_increasing, show.plot = FALSE)
attributes(ItemsList)$intersections$WT









###Increasing in all dust samples
shared_all= significantDG[rownames(significantDG) %in% Increased_shared_all,]
###Increasing only in WT (10 most abundant)
Increased_WT = significantDG[rownames(significantDG) %in% attributes(ItemsList)$intersections$WT,]
Increased_WT =   Increased_WT[order(-Increased_WT$baseMean),] 
Increased_WT = Increased_WT[1:10,]
###Increasing only in 500 (10 most abundant)
Increased_500 = significant500[rownames(significant500) %in% attributes(ItemsList)$intersections$`Coarse SS`,]
Increased_500 =   Increased_500[order(-Increased_500$baseMean),] 
Increased_500 = Increased_500[1:10,]
###Increasing only in 75 (10 most abundant)
Increased_75 = significant75[rownames(significant75) %in% attributes(ItemsList)$intersections$`Medium SS`,]
Increased_75 =   Increased_75[order(-Increased_75$baseMean),] 
Increased_75 = Increased_75[1:10,]
###Increasing only in 25 (10 most abundant)
Increased_25 = significant25[rownames(significant25) %in% attributes(ItemsList)$intersections$`Fine SS`,]
Increased_25 =   Increased_25[order(-Increased_25$baseMean),] 
Increased_25 = Increased_25[1:10,]











####Decreasing######
decreasing500 = rownames(significant500[which(significant500$log2FoldChange < 0 ),]) 
decreasing75  = rownames(significant75[which(significant75$log2FoldChange < 0 ),])
decreasing25  = rownames(significant25[which(significant25$log2FoldChange < 0 ),])
decreasingDG   = rownames(significantDG[which(significantDG$log2FoldChange < 0 ),])

decreased_shared_all = Reduce(intersect, list(decreasing500,decreasing75,decreasing25,decreasingDG ) )#none
decreased_shared_all #[1] "ASVfun21"  "ASVfun41"  "ASVfun411" "ASVfun544" "ASVfun563" "ASVfun768" "ASVfun998"
list_decreasing = list(`Coarse LP` =decreasing500, `Medium LP` = decreasing75, `Fine LP` =  decreasing25, `WT` = decreasingDG )






#####Venn Diagram
b  = ggvenn(
  list_decreasing, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
ggarrange(a,b, labels = "AUTO")











######decreasing 
shared_all= significantDG[rownames(significantDG) %in% decreased_shared_all,]
shared_all =   shared_all[order(-shared_all$baseMean),] 
shared_all = shared_all[1:10,]





