library(FSA)
library(rcompanion)
library(ggpubr)
library(vegan)




####K r strategists
ASVtab_dust<-readRDS("Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
sam_data_dust = data.frame(ASVtab_dust@sam_data) #extract dust data from phyloseq object
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "Dust_Gen"] <- "WT"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "500"] <- "Coarse SS"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "75"] <- "Medium SS"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "25"] <- "Fine SS"
sam_data_dust$Dust_Type<- factor(sam_data_dust$Dust_Type,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))
asvtab_dust = data.frame(ASVtab_dust@otu_table)
asvtab_dust = decostand(asvtab_dust, method = "total", MARGIN = 1)
asvtab_dust = data.frame(t(asvtab_dust))
taxa = data.frame(ASVtab_dust@tax_table)
asvtab_dust$phylum = taxa$Phylum

aggregated = aggregate(. ~phylum, asvtab_dust, FUN = sum)
rownames(aggregated) = aggregated[,1]
aggregated = aggregated[,-1]
aggregated = data.frame(t(aggregated))


Sum_k_strateigists = aggregated$Acidobacteria + aggregated$Actinobacteria
Sum_r_strateigists = aggregated$Proteobacteria + aggregated$Bacteroidetes
ratio1 = Sum_k_strateigists/Sum_r_strateigists

ratio1 = data.frame(ratio = ratio1, sample_type = sam_data_dust$Dust_Type)

kruskal.test(ratio ~ sample_type, data=ratio1)

#######kruskal wallis test
#Kruskal-Wallis chi-squared = 20.148, df = 4, p-value = 0.0004669)
DT = FSA::dunnTest(ratio ~ sample_type,
              data=ratio1,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)
k_r_ratio = ggplot(ratio1, aes(x=sample_type, y=ratio)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","b","b","b","b" ),fun = max,  vjust = 0) +ylab("Ratio") +xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ggtitle("Bact/Arch AA/PB ratio")
k_r_ratio
saveRDS(k_r_ratio, "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4panelE.RDS")







########Fungi
ASVtab_dust<-readRDS("Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
sam_data_dust = data.frame(ASVtab_dust@sam_data) #extract dust data from phyloseq object
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "Dust_Gen"] <- "WT"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "500"] <- "Coarse SS"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "75"] <- "Medium SS"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "25"] <- "Fine SS"
sam_data_dust$Dust_Type<- factor(sam_data_dust$Dust_Type,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))
asvtab_dust = data.frame(ASVtab_dust@otu_table)
asvtab_dust = decostand(asvtab_dust, method = "total", MARGIN = 1)
asvtab_dust = data.frame(t(asvtab_dust))
taxa = data.frame(ASVtab_dust@tax_table)
asvtab_dust$phylum = taxa$Phylum

aggregated = aggregate(. ~phylum, asvtab_dust, FUN = sum)
rownames(aggregated) = aggregated[,1]
aggregated = aggregated[,-1]
aggregated = data.frame(t(aggregated))


Sum_k_strateigists = aggregated$p__Basidiomycota
Sum_r_strateigists = aggregated$p__Ascomycota
ratio1 = Sum_k_strateigists/Sum_r_strateigists

ratio1 = data.frame(ratio = ratio1, sample_type = sam_data_dust$Dust_Type)

kruskal.test(ratio ~ sample_type, data=ratio1)

#######kruskal wallis test
#Kruskal-Wallis chi-squared = 20.148, df = 4, p-value = 0.0004669)
DT = dunnTest(ratio ~ sample_type,
              data=ratio1,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)
k_r_ratio_fung = ggplot(ratio1, aes(x=sample_type, y=ratio)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ),fun = max,  vjust = 0) +ylab("Ratio") +xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ggtitle("Fungi B/A ratio")
k_r_ratio_fung
saveRDS(k_r_ratio_fung, "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4panelH.RDS")


a = ggarrange(k_r_ratio, k_r_ratio_fung, ncol = 1)
a
k_r_ratio
k_r_ratio_fung
ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4/panelD.pdf",   device = pdf)


