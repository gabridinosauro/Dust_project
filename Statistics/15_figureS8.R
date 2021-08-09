#### Load the data for soil###
ASVtab_soil=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soil16sASVnames.RDS")
sam_data_soil = data.frame(ASVtab_soil@sam_data) #extract soil data from phyloseq object
ASVtab_soil <- prune_taxa(taxa_sums(ASVtab_soil) > 0, ASVtab_soil)
ASVtab_soil_merged <- merge_samples(ASVtab_soil, sample_data(ASVtab_soil)$point, fun=mean) ####merge togheter samples at each location
ASVtab_soil_merged = phyloseq_transform_css(ASVtab_soil_merged)
asvtab_soil_merged = data.frame(t(ASVtab_soil_merged@otu_table))  #extract ASV table merged
geodata  = readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_data.RDS")
#####calculate niche breadth
### use the formula from https://doi.org/10.1093/femsec/fiw174
# niche breadth = 1 / sum (average proportion of ASV  per each point ^2) not average abundance!!!!
# calculate relative abundance of each SV per sample
asvtab_soil_dec = decostand(t(asvtab_soil_merged), method = "total")
niche_breadth = asvtab_soil_dec^2
niche_breadth1 = rowSums(niche_breadth)
geodata$nb = 1/niche_breadth1
boxplot(geodata$nb)
plot(geodata$OC_points_percent~geodata$nb)
library(ggplot2)
bac = ggplot( geodata, aes(OC_points_percent,nb)) + geom_point() +
        ylab("Niche breadth") + xlab("Occupancy %") +
        theme_classic() + ggtitle("Bact/Arch") + stat_cor(method = "pearson", label.x = 3, label.y = 30)
bac




ASVtab_soil=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soilITS_ASVnames.RDS")
sam_data_soil = data.frame(ASVtab_soil@sam_data) #extract soil data from phyloseq object
#ASVtab_soil = prune_samples(sample_names(ASVtab_soil)!="DP24-01-S", ASVtab_soil) #eliminate point 24 (low reads)
ASVtab_soil <- prune_taxa(taxa_sums(ASVtab_soil) > 0, ASVtab_soil)
ASVtab_soil_merged <- merge_samples(ASVtab_soil, sample_data(ASVtab_soil)$point, fun=mean) ####merge togheter samples at each location
ASVtab_soil_merged = phyloseq_transform_css(ASVtab_soil_merged)
asvtab_soil_merged = data.frame(t(ASVtab_soil_merged@otu_table))  #extract ASV table merged
geodata  = readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS")
#####calculate niche breadth
### use the formula from https://doi.org/10.1093/femsec/fiw174
# niche breadth = 1 / sum (average proportion of ASV  per each point ^2) not average abundance!!!!
# calculate relative abundance of each SV per sample
asvtab_soil_dec = decostand(t(asvtab_soil_merged), method = "total")
niche_breadth = asvtab_soil_dec^2
niche_breadth1 = rowSums(niche_breadth)
geodata$nb = 1/niche_breadth1
boxplot(geodata$nb)
plot(geodata$OC_points_percent~geodata$nb)
library(ggplot2)
fun = ggplot( geodata, aes(OC_points_percent,nb)) + geom_point() +
  ylab("Niche breadth") + xlab("Occupancy %") +
  theme_classic()+ ggtitle("Fungi") + stat_cor(method = "pearson", label.x = 3, label.y = 30)

library(ggpubr)
ggarrange(bac, fun, labels = "AUTO", nrow = 1)
