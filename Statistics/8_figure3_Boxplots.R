####Aggregate by genus and then run indicator species analyses.
library(DESeq2)
library(vegan)
library(ggplot)
#### 
###Boxplots
ASVtab_dust<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
sam_data_dust = data.frame(ASVtab_dust@sam_data) #extract dust data from phyloseq object
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "Dust_Gen"] <- "WT"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "500"] <- "Coarse LP"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "75"] <- "Medium LP"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "25"] <- "Fine LP"
sam_data_dust$Dust_Type<- factor(sam_data_dust$Dust_Type,levels = c("Soil", "Coarse LP", "Medium LP", "Fine LP", "WT"))
asvtab_dust = data.frame(ASVtab_dust@otu_table) #extract ASVt table not transformed
asvtab_dust = decostand(asvtab_dust, method = "total")
tax_tab = data.frame(ASVtab_dust@tax_table) 

#tax_tab_mas = tax_tab[which(tax_tab$Genus == "Microvirga"),]









######Massilia ASV19
Massilia =  data.frame(asvtab_dust$ASV19)
rownames(asvtab_dust) == rownames(sam_data_dust) #ok
Massilia$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASV19
DEseq2 = readRDS( file = "/Users/gabri/Documents/GitHub/Dust_project/data/Deseq2/deseq_res.RDS")
resultsNames(DEseq2)

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV19",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV19",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV19",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV19",] #SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV19",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV19",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV19",] #NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV19",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV19",] #NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV19",] #NS


outliers = data.frame(sample_type = c("WT"),outlier_txt = c("0.085"))


ASV19_Massilia = ggplot(Massilia, aes(x=sample_type, y=asvtab_dust.ASV19)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","ab","ab","ab","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x = element_blank()) +ggtitle("ASV19, Massilia")  + coord_cartesian(ylim=c(0,0.041)) +
  geom_segment(data = outliers, aes(y=0.041*0.95,yend=0.041,
               xend=factor(sample_type)),
               arrow = arrow(length = unit(0.1,"cm"))) +
  geom_text(data=outliers,aes(y=0.041,label=outlier_txt),
            size=3,vjust=2.8)
  
ASV19_Massilia









###DEseq2_results ### ASV15_Microvirga
microvirga =  data.frame(asvtab_dust$ASV15)
microvirga$sample_type = sam_data_dust$Dust_Type

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV15",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV15",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV15",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV15",]#SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV15",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV15",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV15",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV15",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV15",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV15",]#NS


ASV15_microvirga = ggplot(microvirga, aes(x=sample_type, y=asvtab_dust.ASV15)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","ab","ab","ab","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASV15, Microvirga") +
  theme(axis.title.y=element_blank())
ASV15_microvirga

















###DEseq2_results ### ASV36_Kocuria
Kocuria =  data.frame(asvtab_dust$ASV36)
Kocuria$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASV19
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV36",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV36",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV36",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV36",]#SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV36",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV36",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV36",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV36",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV36",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV36",]#NS


ASV36_Kocuria = ggplot(Kocuria, aes(x=sample_type, y=asvtab_dust.ASV36)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","ab","ab","ab","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASV36, Kocuria") +
  theme(axis.title.y=element_blank())
ASV36_Kocuria












###DEseq2_results ### ASV13_Rubellimicrobium
Rubellimicrobium =  data.frame(asvtab_dust$ASV13)
Rubellimicrobium$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASV13
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV13",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV13",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV13",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV13",]#SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV13",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV13",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV13",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV13",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV13",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV13",]#NS



ASV13_Rubellimicrobium = ggplot(Rubellimicrobium, aes(x=sample_type, y=asvtab_dust.ASV13)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","b","b","b","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASV13, Rubellimicrobium")
ASV13_Rubellimicrobium












###DEseq2_results ### ASV101_Blastococcus
Blastococcus =  data.frame(asvtab_dust$ASV101)
Blastococcus$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASV101
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05))  ; res["ASV101",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05))  ; res["ASV101",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05))  ; res["ASV101",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05))  ; res["ASV101",]#NS



ASV101_Blastococcus = ggplot(Blastococcus, aes(x=sample_type, y=asvtab_dust.ASV101)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","b","b","b","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x = element_blank()) +ggtitle("ASV101, Blastococcus") +
  theme(axis.title.y=element_blank())
ASV101_Blastococcus











###DEseq2_results ### ASV57_Geodermatophilus
Geodermatophilus =  data.frame(asvtab_dust$ASV57)
Geodermatophilus$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASV57
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV57",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV57",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV57",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV57",]#SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) #1.0000000
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) #1.0000000
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) #1.0000000

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) #1.0000000
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05))#1.0000000

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) #7.265133e-01

outliers = data.frame(sample_type = c("Coarse LP", "Medium LP", "Fine LP"), outlier_txt = c("0.026", "0.021", "0.020"))

ASV57_Geodermatophilus = ggplot(Geodermatophilus, aes(x=sample_type, y=asvtab_dust.ASV57)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","b","b","b","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASV57, Geodermatophilus") +  coord_cartesian(ylim=c(0,0.018)) +
  geom_segment(data = outliers, aes(y=0.018*0.95,yend=0.018,
                                    xend=factor(sample_type)),
               arrow = arrow(length = unit(0.1,"cm"))) +
  geom_text(data=outliers,aes(y=0.018,label=outlier_txt),
            size=3,vjust=2.8)+
  theme(axis.title.y=element_blank())
ASV57_Geodermatophilus
















###DEseq2_results ### ASV118_Acidimicrobiia
Acidimicrobiia =  data.frame(asvtab_dust$ASV118)
Acidimicrobiia$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASV101
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05))  ; res["ASV101",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05))  ; res["ASV101",]#SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV101",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV101",]#NS



ASV118_Acidimicrobiia = ggplot(Acidimicrobiia, aes(x=sample_type, y=asvtab_dust.ASV118)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","b","b","b","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASV118, Acidimicrobiia")
ASV118_Acidimicrobiia












###DEseq2_results ### ASV201_Solirubrobacter
Streptomyces =  data.frame(asvtab_dust$ASV201)
Streptomyces$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASV201
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05))  ; res["ASV201",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV201",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV201",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV201",]#SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV201",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV201",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV201",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV201",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV201",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV201",]#NS



ASV201_Streptomyces = ggplot(Streptomyces, aes(x=sample_type, y=asvtab_dust.ASV201)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","ab","ab","ab","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASV201, Streptomyces")+
  theme(axis.title.y=element_blank())
ASV201_Streptomyces
















###DEseq2_results ### ASV300_Solirubrobacter
Solirubrobacter =  data.frame(asvtab_dust$ASV300)
Solirubrobacter$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASV300
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05))  ; res["ASV300",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV300",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV300",]#SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV300",]#SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV300",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV300",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV300",]#NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV300",]#NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASV300",]#SIG

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASV300",]#SIG



ASV300_Solirubrobacter = ggplot(Solirubrobacter, aes(x=sample_type, y=asvtab_dust.ASV300)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","b","b","b","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASV300, Solirubrobacter")+
  theme(axis.title.y=element_blank())
ASV300_Solirubrobacter












####FUNGI#####


#### 
###Boxplots
ASVtab_dust<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
sam_data_dust = data.frame(ASVtab_dust@sam_data) #extract dust data from phyloseq object
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "Dust_Gen"] <- "WT"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "500"] <- "Coarse LP"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "75"] <- "Medium LP"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "25"] <- "Fine LP"
sam_data_dust$Dust_Type<- factor(sam_data_dust$Dust_Type,levels = c("Soil", "Coarse LP", "Medium LP", "Fine LP", "WT"))
asvtab_dust = data.frame(ASVtab_dust@otu_table) #extract ASVt table not transformed
asvtab_dust = decostand(asvtab_dust, method = "total")
tax_tab = data.frame(ASVtab_dust@tax_table) 

#tax_tab_mas = tax_tab[which(tax_tab$Genus == "Microvirga"),]









######Naganishia ASVfun9
Naganishia =  data.frame(asvtab_dust$ASVfun9)
rownames(asvtab_dust) == rownames(sam_data_dust) #ok
Naganishia$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASVfun9
DEseq2 = readRDS( file = "/Users/gabri/Documents/GitHub/Dust_project/Bioinformatics/data/DEseq2_fun.RDS")
resultsNames(DEseq2)

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun9",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun9",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun9",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun9",] #SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun9",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun9",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun9",] #NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun9",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun9",] #NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun9",] #SIG


ASV19_Naganishia = ggplot(Naganishia, aes(x=sample_type, y=asvtab_dust.ASVfun9)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","ab","ab","ab","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASVfun9, Naganishia")
ASV19_Naganishia




######Filobasidium ASVfun13
Filobasidium =  data.frame(asvtab_dust$ASVfun13)
rownames(asvtab_dust) == rownames(sam_data_dust) #ok
Filobasidium$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASVfun9
DEseq2 = readRDS( file = "/Users/gabri/Documents/GitHub/Dust_project/Bioinformatics/data/DEseq2_fun.RDS")
resultsNames(DEseq2)

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun13",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun13",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun13",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun13",] #SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun13",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun13",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun13",] #SIG

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun13",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun13",] #SIG

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun13",] #SIG

outliers = data.frame(sample_type = c("WT", "Fine LP", "Coarse LP"),outlier_txt = c("0.025", "0.10", "0.08"))

ASV13_Filobasidium = ggplot(Filobasidium, aes(x=sample_type, y=asvtab_dust.ASVfun13)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","a","a","a","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASVfun13, Filobasidium") + coord_cartesian(ylim=c(0,0.06)) +
  geom_segment(data = outliers, aes(y=0.06*0.95,yend=0.06,
                                    xend=factor(sample_type)),
               arrow = arrow(length = unit(0.1,"cm"))) +
  geom_text(data=outliers,aes(y=0.06,label=outlier_txt),
            size=3,vjust=2.8) +
  theme(axis.title.y=element_blank())

ASV13_Filobasidium















######Neocamarosporium ASVfun5
Mortierella =  data.frame(asvtab_dust$ASVfun5)
rownames(asvtab_dust) == rownames(sam_data_dust) #ok
Mortierella$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASVfun5
DEseq2 = readRDS( file = "/Users/gabri/Documents/GitHub/Dust_project/Bioinformatics/data/DEseq2_fun.RDS")
resultsNames(DEseq2)

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun5",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun5",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun5",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun5",] #SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun5",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun5",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun5",] #SIG

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun5",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun5",] #SIG

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun5",] #SIG



ASV5_Neocamarosporium = ggplot(Mortierella, aes(x=sample_type, y=asvtab_dust.ASVfun5)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","a","a","a","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASVfun5, Neocamarosporium") +
  theme(axis.title.y=element_blank())
ASV5_Neocamarosporium










######Subramaniula  ASVfun21
Subramaniula  =  data.frame(asvtab_dust$ASVfun21)
rownames(asvtab_dust) == rownames(sam_data_dust) #ok
Subramaniula $sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASVfun21
DEseq2 = readRDS( file = "/Users/gabri/Documents/GitHub/Dust_project/Bioinformatics/data/DEseq2_fun.RDS")
resultsNames(DEseq2)

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun21",] #SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun21",] #SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun21",] #SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun21",] #SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun21",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun21",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun21",] #NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun21",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun21",] #NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun21",]  #NS

outliers = data.frame(sample_type = c("Soil", "Fine LP", "WT"),outlier_txt = c("0.21","0.14", "0.04"))
## Groups = a, ab , ab ,ab, b
ASV21_Subramaniula  = ggplot(Subramaniula , aes(x=sample_type, y=asvtab_dust.ASVfun21)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","ab","ab","ab","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASVfun21, Subramaniula ") + coord_cartesian(ylim=c(0,0.03)) +
  geom_segment(data = outliers, aes(y=0.03*0.95,yend=0.03,xend=factor(sample_type)),arrow = arrow(length = unit(0.1,"cm"))) +
  geom_text(data=outliers,aes(y=0.03,label=outlier_txt),size=3,vjust=2.8) + geom_boxplot(outlier.shape = NA, alpha = 0)+
  theme(axis.title.y=element_blank())
ASV21_Subramaniula 








######Mortierella ASVfun373
Mortierella =  data.frame(asvtab_dust$ASVfun373)
rownames(asvtab_dust) == rownames(sam_data_dust) #ok
Mortierella$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASVfun373
DEseq2 = readRDS( file = "/Users/gabri/Documents/GitHub/Dust_project/Bioinformatics/data/DEseq2_fun.RDS")
resultsNames(DEseq2)

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun373",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun373",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun373",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun373",] #SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun373",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun373",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun373",] #NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun373",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun373",] #NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun373",] #NS

outliers = data.frame(sample_type = c("Soil"),outlier_txt = c("0.01"))


ASV373_Mortierella = ggplot(Mortierella, aes(x=sample_type, y=asvtab_dust.ASVfun373)) +  geom_point(position = "jitter", alpha = 0.3) +  theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","ab","ab","ab","ab" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASVfun373, Mortierella") + coord_cartesian(ylim=c(0,0.007)) +
  geom_segment(data = outliers, aes(y=0.007*0.95,yend=0.007,xend=factor(sample_type)),arrow = arrow(length = unit(0.1,"cm"))) +
  geom_text(data=outliers,aes(y=0.007,label=outlier_txt),size=3,vjust=2.8) + geom_boxplot(outlier.shape = NA, alpha = 0) +
  theme(axis.title.y=element_blank())
ASV373_Mortierella









######Aspergillus ASVfun41
Aspergillus =  data.frame(asvtab_dust$ASVfun41)
rownames(asvtab_dust) == rownames(sam_data_dust) #ok
Aspergillus$sample_type = sam_data_dust$Dust_Type

###DEseq2_results ### ASVfun41
DEseq2 = readRDS( file = "/Users/gabri/Documents/GitHub/Dust_project/Bioinformatics/data/DEseq2_fun.RDS")
resultsNames(DEseq2)

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","500"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun41",] #SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun41",] #SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","25"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun41",] #SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","Soil","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun41",] #SIG


res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","75"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun41",] #SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun41",] #SIG
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","500","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun41",] #SIG

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","25"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun41",] #NS
res = data.frame(results(DEseq2,  contrast=c("Dust_Type","75","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)); res["ASVfun41",] #NS

res = data.frame(results(DEseq2,  contrast=c("Dust_Type","25","Dust_Gen"), pAdjustMethod = "BH", alpha =  0.05)) ; res["ASVfun41",] #NS

## Groups = a, ab , ab ,ab, b
ASVfun41_Aspergillus = ggplot(Aspergillus, aes(x=sample_type, y=asvtab_dust.ASVfun41)) +  geom_point(position = "jitter", alpha = 0.3) +  theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","ab","ab","ab","ab" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x =  element_blank()) +ggtitle("ASVfun41, Aspergillus") + geom_boxplot(outlier.shape = NA, alpha = 0)
ASVfun41_Aspergillus




library(ggpubr)
a = ggarrange(ASV19_Massilia, ASV15_microvirga, ASV36_Kocuria, 
          ASV13_Rubellimicrobium, ASV101_Blastococcus, ASV57_Geodermatophilus, 
          ASV118_Acidimicrobiia, ASV201_Streptomyces, ASV300_Solirubrobacter,
          ASV19_Naganishia, ASV13_Filobasidium, ASV5_Neocamarosporium, 
          ASVfun41_Aspergillus, ASV373_Mortierella, ASV21_Subramaniula, nrow = 5, ncol = 3, align = "v")
#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/Drafts/Review/new_figure3/boxplots.pdf",a , height = 8.5)

