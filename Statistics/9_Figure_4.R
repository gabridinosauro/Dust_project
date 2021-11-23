library(reshape)
library(ggplot2)
library(ggpubr)



######Analysis significant values DESEq2, 500um ##############################
geodata  = readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_data.RDS")
####Bacteria, 500####
significant500 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21500.RDS")
######store all the varaibles in a new dataframe called graph500ab_oc
graph500ab_oc=data.frame(matrix(nrow = nrow(significant500), ncol = 5))
colnames(graph500ab_oc)=c("Log2(FoldChange)","Occupancy","log(Mean Rel Abundance)", "Geographic Range (km2)", "Max Distance (km)")
graph500ab_oc$`Log2(FoldChange)`=significant500$log2FoldChange
rownames(graph500ab_oc)=rownames(significant500)
graph500ab_oc$Occupancy = geodata[rownames(graph500ab_oc), 4] ######occupancy, per point
graph500ab_oc$`log(Mean Rel Abundance)` = log(geodata[rownames(graph500ab_oc), 6]) ######abundance, per point
graph500ab_oc$`Geographic Range (km2)` = geodata[rownames(graph500ab_oc), 9] ######range
graph500ab_oc$`Max Distance (km)` = geodata[rownames(graph500ab_oc), 8] ######distance


####classify each ASV as increasing and decreasing
graph500ab_oc$delta = ifelse(graph500ab_oc$`Log2(FoldChange)`>0 , "Increase", "Decrease")
###melt the dataframe for ggplot
melted_graph500ab_oc=melt( graph500ab_oc )
#####eliminate the log two variable data (not needed in the plots)
melted_graph500ab_oc=melted_graph500ab_oc[which(melted_graph500ab_oc$variable != "Log2(FoldChange)"),]

#######make a partial violin plot
ggplot(melted_graph500ab_oc, aes(x=delta, y=value)) + 
  geom_violin(trim=FALSE) + facet_wrap(variable ~ ., scales = "free_y") +
  geom_jitter(shape=16, position=position_jitter(0.2))




######bacteria75#####
significant75 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2175.RDS")
######store all the varaibles in a new dataframe called graph75ab_oc
graph75ab_oc=data.frame(matrix(nrow = nrow(significant75), ncol = 5))
colnames(graph75ab_oc)=c("Log2(FoldChange)","Occupancy","log(Mean Rel Abundance)", "Geographic Range (km2)", "Max Distance (km)")
graph75ab_oc$`Log2(FoldChange)`=significant75$log2FoldChange
rownames(graph75ab_oc)=rownames(significant75)
graph75ab_oc$Occupancy = geodata[rownames(graph75ab_oc), 4] ######occupancy, per point
graph75ab_oc$`log(Mean Rel Abundance)` = log(geodata[rownames(graph75ab_oc), 6]) ######abundance, per point
graph75ab_oc$`Geographic Range (km2)` = geodata[rownames(graph75ab_oc), 9] ######range
graph75ab_oc$`Max Distance (km)` = geodata[rownames(graph75ab_oc), 8] ######distance


####classify each ASV as increasing and decreasing
graph75ab_oc$delta = ifelse(graph75ab_oc$`Log2(FoldChange)`>0 , "Increase", "Decrease")
###melt the dataframe for ggplot
melted_graph75ab_oc=melt( graph75ab_oc )
#####eliminate the log two variable data (not needed in the plots)
melted_graph75ab_oc=melted_graph75ab_oc[which(melted_graph75ab_oc$variable != "Log2(FoldChange)"),]

#######make a partial violin plot
ggplot(melted_graph75ab_oc, aes(x=delta, y=value)) + 
  geom_violin(trim=FALSE) + facet_wrap(variable ~ ., scales = "free_y") +
  geom_jitter(shape=16, position=position_jitter(0.2))


######bacteria25#####
significant25 = readRDS(file = "//Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2125.RDS")
######store all the varaibles in a new dataframe called graph25ab_oc
graph25ab_oc=data.frame(matrix(nrow = nrow(significant25), ncol = 5))
colnames(graph25ab_oc)=c("Log2(FoldChange)","Occupancy","log(Mean Rel Abundance)", "Geographic Range (km2)", "Max Distance (km)")
graph25ab_oc$`Log2(FoldChange)`=significant25$log2FoldChange
rownames(graph25ab_oc)=rownames(significant25)
graph25ab_oc$Occupancy = geodata[rownames(graph25ab_oc), 4] ######occupancy, per point
graph25ab_oc$`log(Mean Rel Abundance)` = log(geodata[rownames(graph25ab_oc), 6])######abundance, per point
graph25ab_oc$`Geographic Range (km2)` = geodata[rownames(graph25ab_oc), 9] ######range
graph25ab_oc$`Max Distance (km)` = geodata[rownames(graph25ab_oc), 8] ######distance


####classify each ASV as increasing and decreasing
graph25ab_oc$delta = ifelse(graph25ab_oc$`Log2(FoldChange)`>0 , "Increase", "Decrease")
###melt the dataframe for ggplot
melted_graph25ab_oc=melt( graph25ab_oc )
#####eliminate the log two variable data (not needed in the plots)
melted_graph25ab_oc=melted_graph25ab_oc[which(melted_graph25ab_oc$variable != "Log2(FoldChange)"),]

#######make a partial violin plot
ggplot(melted_graph25ab_oc, aes(x=delta, y=value)) + 
  geom_violin(trim=FALSE) + facet_wrap(variable ~ ., scales = "free_y") +
  geom_jitter(shape=16, position=position_jitter(0.2))





######DG#######
significantDG = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21DG.RDS")
######store all the varaibles in a new dataframe called graphDGab_oc
graphDGab_oc=data.frame(matrix(nrow = nrow(significantDG), ncol = 5))
colnames(graphDGab_oc)=c("Log2(FoldChange)","Occupancy","log(Mean Rel Abundance)", "Geographic Range (km2)", "Max Distance (km)")
graphDGab_oc$`Log2(FoldChange)`=significantDG$log2FoldChange
rownames(graphDGab_oc)=rownames(significantDG)
graphDGab_oc$Occupancy = geodata[rownames(graphDGab_oc), 4] ######occupancy, per point
graphDGab_oc$`log(Mean Rel Abundance)` = log(geodata[rownames(graphDGab_oc), 6]) ######abundance, per point
graphDGab_oc$`Geographic Range (km2)` = geodata[rownames(graphDGab_oc), 9] ######range
graphDGab_oc$`Max Distance (km)` = geodata[rownames(graphDGab_oc), 8] ######distance


####classify each ASV as increasing and decreasing
graphDGab_oc$delta = ifelse(graphDGab_oc$`Log2(FoldChange)`>0 , "Increase", "Decrease")
###melt the dataframe for ggplot
melted_graphDGab_oc=melt( graphDGab_oc )
#####eliminate the log two variable data (not needed in the plots)
melted_graphDGab_oc=melted_graphDGab_oc[which(melted_graphDGab_oc$variable != "Log2(FoldChange)"),]

#######make a partial violin plot
ggplot(melted_graphDGab_oc, aes(x=delta, y=value)) + 
  geom_violin(trim=FALSE) + facet_wrap(variable ~ ., scales = "free_y") +
  geom_jitter(shape=16, position=position_jitter(0.2))



#######plot altogheter#####
melted_graph500ab_oc$Dust_fraction="Coarse SS" 
melted_graph75ab_oc$Dust_fraction="Medium SS" 
melted_graphDGab_oc$Dust_fraction="WT" 
melted_graph25ab_oc$Dust_fraction="Fine SS" 

graph_total = rbind (melted_graph500ab_oc,melted_graph75ab_oc,melted_graphDGab_oc, melted_graph25ab_oc)
graph_total$Dust_fraction <- factor(graph_total$Dust_fraction , levels=c("Coarse SS", "Medium SS", "Fine SS", "WT"))


give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}




graphtotal = ggplot(graph_total, aes(x=Dust_fraction, y=value, fill = delta )) + 
  geom_point(position=position_jitterdodge(),alpha=0.2, aes(fill=delta)) +
  geom_boxplot(trim=FALSE, alpha = 0.6, outlier.shape = NA) + facet_wrap(variable ~ ., scales = "free_y", nrow = 1) +
  theme_classic()  +
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","#4DBBD5B2"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab(NULL) + xlab(NULL) +ggtitle("Bacteria/Archaea")
graphtotal

#+
#  stat_summary(fun.data = give.n, geom = "text", fun.y = median,
#               position = position_dodge(width = 0.75))
##########statistical tests##########
wilcox.test(graph500ab_oc$Occupancy ~  graph500ab_oc$delta)#W = 1040, p-value = 0.2666
wilcox.test(graph75ab_oc$Occupancy ~  graph75ab_oc$delta)#W = 1554.5, p-value = 0.9418
wilcox.test(graph25ab_oc$Occupancy ~  graph25ab_oc$delta)#W = 2067.5, p-value = 0.4531
wilcox.test(graphDGab_oc$Occupancy ~  graphDGab_oc$delta)#W = 3287, p-value = 0.02394

wilcox.test(graph500ab_oc$`log(Mean Rel Abundance)` ~  graph500ab_oc$delta)#W = 1194, p-value = 0.9859
wilcox.test(graph75ab_oc$`log(Mean Rel Abundance)`  ~  graph75ab_oc$delta)#W = 1768, p-value = 0.2444
wilcox.test(graph25ab_oc$`log(Mean Rel Abundance)`  ~  graph25ab_oc$delta)#W = 2613, p-value = 0.0005275
wilcox.test(graphDGab_oc$`log(Mean Rel Abundance)`  ~  graphDGab_oc$delta)#W = 3741, p-value = 0.3237

wilcox.test(graph500ab_oc$`Geographic Range (km2)` ~  graph500ab_oc$delta)#W = 336, p-value = 0.03368
wilcox.test(graph75ab_oc$`Geographic Range (km2)`   ~  graph75ab_oc$delta)#W = 484, p-value = 0.1883
wilcox.test(graph25ab_oc$`Geographic Range (km2)`   ~  graph25ab_oc$delta)#W = 452, p-value = 0.6417
wilcox.test(graphDGab_oc$`Geographic Range (km2)`   ~  graphDGab_oc$delta)#W = 1049.5, p-value = 0.008294


wilcox.test(graph500ab_oc$`Max Distance (km)` ~  graph500ab_oc$delta)#W = 526, p-value = 0.006226
wilcox.test(graph75ab_oc$`Max Distance (km)`   ~  graph75ab_oc$delta)#W = 762, p-value = 0.01446
wilcox.test(graph25ab_oc$`Max Distance (km)`   ~  graph25ab_oc$delta)#W = 875.5, p-value = 0.06603
wilcox.test(graphDGab_oc$`Max Distance (km)`   ~  graphDGab_oc$delta)#W = 2140.5, p-value = 0.003595

library(reshape)
library(ggplot2)





unique(graph_total$variable)
#######eliminate relative abundance
graph_total = graph_total[graph_total$variable != "log(Mean Rel Abundance)",]
graphtotal = ggplot(graph_total, aes(x=Dust_fraction, y=value, fill = delta )) + 
  #geom_point(position=position_jitterdodge(),alpha=0.2, aes(fill=delta)) +
  geom_boxplot(trim=FALSE, alpha = 0.6, outlier.shape = NA) + facet_wrap(variable ~ ., scales = "free_y", nrow = 1) +
  theme_classic()  +
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","#4DBBD5B2"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab(NULL) + xlab(NULL) +ggtitle("Bacteria/Archaea")
graphtotal













######Fungi##############################
geodata  = readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS")
####Bacteria, 500####
significant500 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21500ITS.RDS")
######store all the varaibles in a new dataframe called graph500ab_oc
graph500ab_oc=data.frame(matrix(nrow = nrow(significant500), ncol = 5))
colnames(graph500ab_oc)=c("Log2(FoldChange)","Occupancy","log(Mean Rel Abundance)", "Geographic Range (km2)", "Max Distance (km)")
graph500ab_oc$`Log2(FoldChange)`=significant500$log2FoldChange
rownames(graph500ab_oc)=rownames(significant500)
graph500ab_oc$Occupancy = geodata[rownames(graph500ab_oc), 4] ######occupancy, per point
graph500ab_oc$`log(Mean Rel Abundance)` = log(geodata[rownames(graph500ab_oc), 6]) ######abundance, per point
graph500ab_oc$`Geographic Range (km2)` = geodata[rownames(graph500ab_oc), 9] ######range
graph500ab_oc$`Max Distance (km)` = geodata[rownames(graph500ab_oc), 8] ######distance


####classify each ASV as increasing and decreasing
graph500ab_oc$delta = ifelse(graph500ab_oc$`Log2(FoldChange)`>0 , "Increase", "Decrease")
###melt the dataframe for ggplot
melted_graph500ab_oc=melt( graph500ab_oc )
#####eliminate the log two variable data (not needed in the plots)
melted_graph500ab_oc=melted_graph500ab_oc[which(melted_graph500ab_oc$variable != "Log2(FoldChange)"),]

#######make a partial violin plot
ggplot(melted_graph500ab_oc, aes(x=delta, y=value)) + 
  geom_violin(trim=FALSE) + facet_wrap(variable ~ ., scales = "free_y") +
  geom_jitter(shape=16, position=position_jitter(0.2))




######bacteria75#####
significant75 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2175_ITS.RDS")
######store all the varaibles in a new dataframe called graph75ab_oc
graph75ab_oc=data.frame(matrix(nrow = nrow(significant75), ncol = 5))
colnames(graph75ab_oc)=c("Log2(FoldChange)","Occupancy","log(Mean Rel Abundance)", "Geographic Range (km2)", "Max Distance (km)")
graph75ab_oc$`Log2(FoldChange)`=significant75$log2FoldChange
rownames(graph75ab_oc)=rownames(significant75)
graph75ab_oc$Occupancy = geodata[rownames(graph75ab_oc), 4] ######occupancy, per point
graph75ab_oc$`log(Mean Rel Abundance)` = log(geodata[rownames(graph75ab_oc), 6]) ######abundance, per point
graph75ab_oc$`Geographic Range (km2)` = geodata[rownames(graph75ab_oc), 9] ######range
graph75ab_oc$`Max Distance (km)` = geodata[rownames(graph75ab_oc), 8] ######distance


####classify each ASV as increasing and decreasing
graph75ab_oc$delta = ifelse(graph75ab_oc$`Log2(FoldChange)`>0 , "Increase", "Decrease")
###melt the dataframe for ggplot
melted_graph75ab_oc=melt( graph75ab_oc )
#####eliminate the log two variable data (not needed in the plots)
melted_graph75ab_oc=melted_graph75ab_oc[which(melted_graph75ab_oc$variable != "Log2(FoldChange)"),]

#######make a partial violin plot
ggplot(melted_graph75ab_oc, aes(x=delta, y=value)) + 
  geom_violin(trim=FALSE) + facet_wrap(variable ~ ., scales = "free_y") +
  geom_jitter(shape=16, position=position_jitter(0.2))


######bacteria25#####
significant25 = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq2125_ITS.RDS")
######store all the varaibles in a new dataframe called graph25ab_oc
graph25ab_oc=data.frame(matrix(nrow = nrow(significant25), ncol = 5))
colnames(graph25ab_oc)=c("Log2(FoldChange)","Occupancy","log(Mean Rel Abundance)", "Geographic Range (km2)", "Max Distance (km)")
graph25ab_oc$`Log2(FoldChange)`=significant25$log2FoldChange
rownames(graph25ab_oc)=rownames(significant25)
graph25ab_oc$Occupancy = geodata[rownames(graph25ab_oc), 4] ######occupancy, per point
graph25ab_oc$`log(Mean Rel Abundance)` = log(geodata[rownames(graph25ab_oc), 6]) ######abundance, per point
graph25ab_oc$`Geographic Range (km2)` = geodata[rownames(graph25ab_oc), 9] ######range
graph25ab_oc$`Max Distance (km)` = geodata[rownames(graph25ab_oc), 8] ######distance


####classify each ASV as increasing and decreasing
graph25ab_oc$delta = ifelse(graph25ab_oc$`Log2(FoldChange)`>0 , "Increase", "Decrease")
###melt the dataframe for ggplot
melted_graph25ab_oc=melt( graph25ab_oc )
#####eliminate the log two variable data (not needed in the plots)
melted_graph25ab_oc=melted_graph25ab_oc[which(melted_graph25ab_oc$variable != "Log2(FoldChange)"),]

#######make a partial violin plot
ggplot(melted_graph25ab_oc, aes(x=delta, y=value)) + 
  geom_violin(trim=FALSE) + facet_wrap(variable ~ ., scales = "free_y") +
  geom_jitter(shape=16, position=position_jitter(0.2))





######DG#######
significantDG = readRDS(file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/sigtabDEseq21DGITS.RDS")
######store all the varaibles in a new dataframe called graphDGab_oc
graphDGab_oc=data.frame(matrix(nrow = nrow(significantDG), ncol = 5))
colnames(graphDGab_oc)=c("Log2(FoldChange)","Occupancy","log(Mean Rel Abundance)", "Geographic Range (km2)", "Max Distance (km)")
graphDGab_oc$`Log2(FoldChange)`=significantDG$log2FoldChange
rownames(graphDGab_oc)=rownames(significantDG)
graphDGab_oc$Occupancy = geodata[rownames(graphDGab_oc), 4] ######occupancy, per point
graphDGab_oc$`log(Mean Rel Abundance)` = log(geodata[rownames(graphDGab_oc), 6]) ######abundance, per point
graphDGab_oc$`Geographic Range (km2)` = geodata[rownames(graphDGab_oc), 9] ######range
graphDGab_oc$`Max Distance (km)` = geodata[rownames(graphDGab_oc), 8] ######distance


####classify each ASV as increasing and decreasing
graphDGab_oc$delta = ifelse(graphDGab_oc$`Log2(FoldChange)`>0 , "Increase", "Decrease")
###melt the dataframe for ggplot
melted_graphDGab_oc=melt( graphDGab_oc )
#####eliminate the log two variable data (not needed in the plots)
melted_graphDGab_oc=melted_graphDGab_oc[which(melted_graphDGab_oc$variable != "Log2(FoldChange)"),]

#######make a partial violin plot
ggplot(melted_graphDGab_oc, aes(x=delta, y=value)) + 
  geom_violin(trim=FALSE) + facet_wrap(variable ~ ., scales = "free_y") +
  geom_jitter(shape=16, position=position_jitter(0.2))



#######plot altogheter#####
melted_graph500ab_oc$Dust_fraction="Coarse SS" 
melted_graph75ab_oc$Dust_fraction="Medium SS" 
melted_graphDGab_oc$Dust_fraction="WT" 
melted_graph25ab_oc$Dust_fraction="Fine SS" 

graph_total = rbind (melted_graph500ab_oc,melted_graph75ab_oc,melted_graphDGab_oc, melted_graph25ab_oc)
graph_total$Dust_fraction <- factor(graph_total$Dust_fraction , levels=c("Coarse SS", "Medium SS", "Fine SS", "WT"))


give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}




graphtotal_fung = ggplot(graph_total, aes(x=Dust_fraction, y=value, fill = delta )) + 
  geom_point(position=position_jitterdodge(),alpha=0.2, aes(fill=delta)) +
  geom_boxplot(trim=FALSE, alpha = 0.6, outlier.shape = NA) + facet_wrap(variable ~ ., scales = "free_y", nrow = 1) +
  theme_classic()  +
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","#4DBBD5B2"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab(NULL) + xlab(NULL) +ggtitle("Fungi")
graphtotal_fung

##########statistical tests##########
wilcox.test(graph500ab_oc$Occupancy ~  graph500ab_oc$delta)#NS
wilcox.test(graph75ab_oc$Occupancy ~  graph75ab_oc$delta)#NS
wilcox.test(graph25ab_oc$Occupancy ~  graph25ab_oc$delta)#NS
wilcox.test(graphDGab_oc$Occupancy ~  graphDGab_oc$delta)#NS

wilcox.test(graph500ab_oc$`log(Mean Rel Abundance)` ~  graph500ab_oc$delta)#W = 10, p-value = 0.9091
wilcox.test(graph75ab_oc$`log(Mean Rel Abundance)`  ~  graph75ab_oc$delta)#W = 50, p-value = 0.0006
wilcox.test(graph25ab_oc$`log(Mean Rel Abundance)`  ~  graph25ab_oc$delta)#W = 64, p-value = 0.1087
wilcox.test(graphDGab_oc$`log(Mean Rel Abundance)`  ~  graphDGab_oc$delta)#W = 258, p-value = 0.4284

wilcox.test(graph500ab_oc$`Geographic Range (km2)` ~  graph500ab_oc$delta)#NS
wilcox.test(graph75ab_oc$`Geographic Range (km2)`   ~  graph75ab_oc$delta)#NS
wilcox.test(graph25ab_oc$`Geographic Range (km2)`   ~  graph25ab_oc$delta)#NS
wilcox.test(graphDGab_oc$`Geographic Range (km2)`   ~  graphDGab_oc$delta)#NS


wilcox.test(graph500ab_oc$`Max Distance (km)` ~  graph500ab_oc$delta)#NS
wilcox.test(graph75ab_oc$`Max Distance (km)`   ~  graph75ab_oc$delta)#NS
wilcox.test(graph25ab_oc$`Max Distance (km)`   ~  graph25ab_oc$delta)#NS
wilcox.test(graphDGab_oc$`Max Distance (km)`   ~  graphDGab_oc$delta)#WNS

graph_total = graph_total[graph_total$variable != "log(Mean Rel Abundance)",]
graphtotal_fung = ggplot(graph_total, aes(x=Dust_fraction, y=value, fill = delta )) + 
  #geom_point(position=position_jitterdodge(),alpha=0.2, aes(fill=delta)) +
  geom_boxplot(trim=FALSE, alpha = 0.6, outlier.shape = NA) + facet_wrap(variable ~ ., scales = "free_y", nrow = 1) +
  theme_classic()  +
  scale_fill_manual(name = "Sample", labels = c("Soil", "Dust"), values = c("#996600","#4DBBD5B2"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab(NULL) + xlab(NULL) +ggtitle("Bacteria/Archaea")
graphtotal_fung




ggarrange(graphtotal,graphtotal_fung, ncol = 1, common.legend = TRUE)
ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/4th_draft/new_figures/figure4.pdf", device = "pdf", height = 6, width = 7.5)





