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
geodata  = readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_data.RDS") # to produce this file run the code "geographic_calculations.RDS"
ASVtab_dust<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
ASVtab_dust=prune_samples(sample_data(ASVtab_dust)$point!=26, ASVtab_dust)
ASVtab_dust <- prune_taxa(taxa_sums(ASVtab_dust) > 0, ASVtab_dust)
sam_data_dust = data.frame(ASVtab_dust@sam_data) #extract dust data from phyloseq object
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "Dust_Gen"] <- "WT"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "500"] <- "Coarse LP"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "75"] <- "Medium LP"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "25"] <- "Fine LP"
sam_data_dust$Dust_Type<- factor(sam_data_dust$Dust_Type,levels = c("Soil", "Coarse LP", "Medium LP", "Fine LP", "WT"))




asvtab_dust = data.frame(ASVtab_dust@otu_table) #extract ASVt table not transformed






####set a threshold for the first analysis
threshold = 2






###select the ASVs based on high occupancy
selection = round((nrow(geodata)/100)*threshold) ###define the number of ASVs to pick
high_occupancy_ASVs = geodata[order(-geodata$OC_points_percent),] ### order the table
high_occupancy_ASVs = high_occupancy_ASVs[1:selection,] ### extract those with the highest occupancy
min = min(high_occupancy_ASVs$OC_points_percent) ####correct for those having the same occupancy
high_occupancy_ASVs = geodata[which(geodata$OC_points_percent >= min),] ##finally select them








###select the ASVs based on high abundance
`%nin%` = Negate(`%in%`)
geodata_specialists = geodata[rownames(geodata) %nin% rownames(high_occupancy_ASVs), ]
selection = round((nrow(geodata_specialists)/100)*threshold) ###define the number of ASVs to pick
high_abundance_ASVs = geodata_specialists[order(-geodata_specialists$Ab_bypoint),] ### order the table
high_abundance_ASVs = high_abundance_ASVs[1:selection,] ### extract those with the highest occupancy
min_ab = min(high_abundance_ASVs$Ab_bypoint) ####correct for those having the same occupancy
high_abundance_ASVs = geodata_specialists[which(geodata_specialists$Ab_bypoint >= min_ab),] ##finally select them





















##########first method############
##Plot  the occupancy- ab curve with colors curve
oc_ab_plot = geodata[ ,c(4,6) ] 
oc_ab_plot$color = NA 
oc_ab_plot[which(rownames(oc_ab_plot) %in% rownames(high_occupancy_ASVs)), 3] = "Generalists"
oc_ab_plot[which(rownames(oc_ab_plot) %in% rownames(high_abundance_ASVs)), 3] = "Specialists"
oc_ab_plot[which(rownames(oc_ab_plot) %nin% rownames(high_occupancy_ASVs) & rownames(oc_ab_plot) %nin% rownames(high_abundance_ASVs)), 3] = "Unclassifed"


ab_oc_curve = ggplot(oc_ab_plot, aes(Ab_bypoint, OC_points_percent)) + 
  geom_point(aes (color= color), size = 1) +
  scale_x_log10() + 
  theme_classic() + 
  xlab("Log(average abundance)") +
  ylab("Occupancy %") +
  theme(legend.title=element_blank(), legend.position = c(0.23, 0.85),
        legend.background = element_rect(fill = "white", color = "black")) +
  ggtitle("Bacteria/Archaea, 2% threshold") +
  scale_color_manual( labels = c("Generalists", "Specialists", "Unclassified"),values=c("#56B4E9", "#E69F00", "#999999")) #+
ab_oc_curve

#pdf("/Users/gabri/OneDrive - University of Arizona/dust_project/3rd_draft/plots/new_figures/ab_oc_curve_bac.pdf", width = 3, height = 3)
ab_oc_curve
#dev.off()
#saveRDS(ab_oc_curve,  "/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/ab_oc_curve_plot_bac.RDS" )

#### plot the abundance occupancy curve
#ab_oc_curve = geodata[,c(4,6)]
#ab_oc_curve_plot_fung = ggplot(ab_oc_curve, aes(Ab_bypoint, OC_points_percent)) + 
#  geom_point(color = "grey")+
#  scale_x_log10() + 
#  geom_hline(yintercept = min -2) +
#  geom_vline(xintercept = min_ab )+
#  theme_classic() + 
#  xlab("Log(average abundance)") +
#  ylab("Occupancy %") +
#  scale_color_manual(values = c( "#E69F00","#999999")) +
#  theme(legend.title=element_blank(), legend.position = c(0.15, 0.85),
#        legend.background = element_rect(fill = "white", color = "black")) +
#  ggtitle("Fungi, 2% thresholds") 
#ab_oc_curve_plot_fung
#saveRDS(ab_oc_curve_plot_bac,  "/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/ab_oc_curve_plot_bac.RDS" )


####calculate cumulative relative abundance
rar_ab_per_point = data.frame(matrix(ncol = 1, nrow = nrow(asvtab_dust)))
colnames(rar_ab_per_point)  = c("high")
rownames(rar_ab_per_point) = rownames(asvtab_dust) 


for (i in 1:nrow(asvtab_dust)) 
{
  sample = t(data.frame(asvtab_dust[i,]))
  sample = data.frame(sample[apply(sample, 1, function(row) all(row !=0 )), ])
  rar_ab_per_point[i,1] = sum(sample[rownames(sample) %in% rownames(high_occupancy_ASVs),]) / sum(sample)
  
}
rar_ab_per_point$Dust_fraction = sam_data_dust$Dust_Type


kruskal.test(high ~ Dust_fraction, data=rar_ab_per_point)
#######kruskal wallis test
#Kruskal-Wallis chi-squared = 12.8, df = 4, p-value = 0.01)
DT = dunnTest(high ~ Dust_fraction,
              data=rar_ab_per_point,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

####plot it with ggplot
melted = melt(rar_ab_per_point, stringsAsFactors = FALSE)

####ggplot
relative_abundance_high_occupancy_bac = ggplot(melted, aes(x=Dust_fraction, y=value)) +  geom_point(position = "jitter",  alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","b","ab","ab","b" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ggtitle("Bact/Arch, Generalists")
relative_abundance_high_occupancy_bac
#saveRDS(relative_abundance_high_occupancy_bac, "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4panelA.RDS")
#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4/panelA.pdf",   device = pdf)



####caclulate high abundance
####calculate cumulative relative abundance
rar_ab_per_point = data.frame(matrix(ncol = 1, nrow = nrow(asvtab_dust)))
colnames(rar_ab_per_point)  = c("high")
rownames(rar_ab_per_point) = rownames(asvtab_dust) 


for (i in 1:nrow(asvtab_dust)) 
{
  sample = t(data.frame(asvtab_dust[i,]))
  sample = data.frame(sample[apply(sample, 1, function(row) all(row !=0 )), ])
  rar_ab_per_point[i,1] = sum(sample[rownames(sample) %in% rownames(high_abundance_ASVs),]) / sum(sample)
  
}
rar_ab_per_point$Dust_fraction = sam_data_dust$Dust_Type



#######kruskal wallis test
#Kruskal-Wallis chi-squared = 20.148, df = 4, p-value = 0.0004669)
DT = dunnTest(high ~ Dust_fraction,
              data=rar_ab_per_point,
              method="bh")      # Adjusts p-values for multiple comparisons;
PT = DT$res
PT
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

####plot it with ggplot
melted = melt(rar_ab_per_point, stringsAsFactors = FALSE)



####ggplot
relative_abundance_high_abudance_bac = ggplot(melted, aes(x=Dust_fraction, y=value)) +  geom_point(position = "jitter", alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ),fun = max,  vjust = 0) +ylab(NULL) +xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ggtitle("Bact/Arch, Specialists") +ylab("Relative abundance")
relative_abundance_high_abudance_bac
#saveRDS(relative_abundance_high_abudance_bac, "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4panelB.RDS")
#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4/panelB.pdf",   device = pdf)

a =ggarrange(ab_oc_curve, relative_abundance_high_occupancy_bac, relative_abundance_high_abudance_bac, ncol = 3)
a
ggsave("gen_spec.pdf", plot = a ,"/Users/gabri/OneDrive - University of Arizona/dust_project/4th_draft/new_figures/", device = "pdf", height = 3 , width = 4.5)
#saveRDS(a, file = "/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/boxplots_bacteria.rds" )



#list = list(ab_oc_curve_plot_bac, relative_abundance_high_occupancy_bac,relative_abundance_high_abudance_bac)

#grid.arrange(
#  grobs = list,
#  layout_matrix = rbind(c(1, 1, 2),
#                        c(1, 1, 3))
#)