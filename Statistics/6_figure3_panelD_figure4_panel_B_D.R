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
geodata  = readRDS("Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS") # to produce this file run the code "geographic_calculations.RDS"
ASVtab_dust<-readRDS("Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
sam_data_dust = data.frame(ASVtab_dust@sam_data) #extract dust data from phyloseq object
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "Dust_Gen"] <- "WT"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "500"] <- "Coarse SP"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "75"] <- "Medium SP"
sam_data_dust$Dust_Type[sam_data_dust$Dust_Type == "25"] <- "Fine SP"
sam_data_dust$Dust_Type<- factor(sam_data_dust$Dust_Type,levels = c("Soil", "Coarse SP", "Medium SP", "Fine SP", "WT"))




asvtab_dust = data.frame(ASVtab_dust@otu_table) #extract ASVt table not transformed






####set a threshold for the first analysis
threshold = 1






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
oc_ab_plot[which(rownames(oc_ab_plot) %in% rownames(high_occupancy_ASVs)), 3] = "High_Occupancy"
oc_ab_plot[which(rownames(oc_ab_plot) %in% rownames(high_abundance_ASVs)), 3] = "High_Abundance"
oc_ab_plot[which(rownames(oc_ab_plot) %nin% rownames(high_occupancy_ASVs) & rownames(oc_ab_plot) %nin% rownames(high_abundance_ASVs)), 3] = "Unclassifed"


ab_oc_curve = ggplot(oc_ab_plot, aes(Ab_bypoint, OC_points_percent)) + 
  geom_point(aes (color= color)) +
  scale_x_log10() + 
  theme_classic() + 
  xlab("Log(average abundance)") +
  ylab("Occupancy %") +
  theme(legend.title=element_blank(), legend.position = c(0.23, 0.85),
        legend.background = element_rect(fill = "white", color = "black")) +
  ggtitle("Fungi, 1% threshold") +
  scale_color_manual( labels = c("Specialists", "Generalists", "Unclassified"),values=c("#56B4E9", "#E69F00", "#999999")) #+
ab_oc_curve

#ggsave( "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/ab_oc_curve_fung.pdf",
#        device = "pdf", 
#        ab_oc_curve)



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
relative_abundance_high_occupancy_bac = ggplot(melted, aes(x=Dust_fraction, y=value)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ),fun = max,  vjust = 0) +ylab("Relative abundance") +xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ggtitle("Fungi, generalists")
relative_abundance_high_occupancy_bac
#saveRDS(relative_abundance_high_occupancy_bac, "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4panelC.RDS")
#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4/panel_11.pdf",   device = pdf)

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
relative_abundance_high_abudance_bac = ggplot(melted, aes(x=Dust_fraction, y=value)) +  geom_point(position = "jitter") +
  geom_boxplot(outlier.shape = NA, alpha = 0)  + theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ),fun = max,  vjust = 0) +ylab(NULL) +xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ggtitle("Fungi, specialists") +ylab("Relative abundance")
relative_abundance_high_abudance_bac
#saveRDS(relative_abundance_high_abudance_bac, "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4panelD.RDS")
#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/figure4/panel_12.pdf",   device = pdf)

a =ggarrange(relative_abundance_high_occupancy_bac, relative_abundance_high_abudance_bac, ncol = 2)
a
saveRDS(a, file = "/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/boxplots_fungi.rds" )



list = list(ab_oc_curve_plot_bac, relative_abundance_high_occupancy_bac,relative_abundance_high_abudance_bac)

grid.arrange(
  grobs = list,
  layout_matrix = rbind(c(1, 1, 2),
                        c(1, 1, 3))
)

