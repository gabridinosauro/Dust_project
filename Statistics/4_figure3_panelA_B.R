####Distance decay Pattern####

library(vegan)#ecology
library(sp)#spatial data
library(metagMisc)#CSS transformation
library(ade4)


###Load the data####
asvtab <- readRDS("Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soil16sASVnames.RDS")
ASVtab_merged= merge_samples(asvtab, sample_data(asvtab)$point, fun=mean) # merge points togheter
asvtab_m <- prune_taxa(taxa_sums(ASVtab_merged) > 0, ASVtab_merged)
asvtab_m = phyloseq_transform_css(asvtab_m)
sam_data = data.frame(asvtab_m@sam_data)
asvtab_m = data.frame(asvtab_m@otu_table)
sam_data = sam_data[,-c(1,2)]



####Load coordinates and calculate a distance matrix
coordinate <- read.csv("Documents/GitHub/Dust_project/data/metadata/points_coordinate.csv")
coordinate=coordinate[order(coordinate$id),]
rownames(coordinate)=coordinate$id
coordinate=coordinate[-26,]
coordinates(coordinate) <- cbind(coordinate$X , coordinate$Y) #### Transform it into a spatial file
proj4string(coordinate) = CRS("+proj=longlat +datum=WGS84 +no_defs")
distmatrix<-geosphere::distm(coordinate) 

beta_div = vegdist(t(asvtab_m), method="bray")
micro.dist = as.matrix(vegdist(t(asvtab_m), method="bray"))
databac <- data.frame(Dissimilarity = as.vector(micro.dist),
                      Distance = as.vector(distmatrix))

databac = databac[which (databac$Dissimilarity!=0),]
cor.test(databac$Dissimilarity, databac$Distance)#nothing



#mantel R
beta_div = vegdist(t(asvtab_m), method="bray")
dist.mat <- as.dist(distmatrix)

mantel.rtest(beta_div, dist.mat, nrepet = 9999)




a = ggplot(databac, aes(Distance,Dissimilarity))+
  geom_point(shape = 1) +
  theme_bw() +
  ggtitle("Bact/Arc") +
  theme_classic()+
  xlab("Distance (m)")

a

#ggsave( "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/dist_dec_bac.pdf",
#        device = "pdf", 
#        a)


#saveRDS(a, "/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/cor_bac.RDS")















###########################FUNGI##################################################
##Load the data####
asvtab <- readRDS("Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soilITS_ASVnames.RDS")
ASVtab_merged= merge_samples(asvtab, sample_data(asvtab)$point, fun=mean) # merge points togheter
asvtab_m <- prune_taxa(taxa_sums(ASVtab_merged) > 0, ASVtab_merged)
asvtab_m = phyloseq_transform_css(asvtab_m)
sam_data = data.frame(asvtab_m@sam_data)
asvtab_m = data.frame(asvtab_m@otu_table)
sam_data = sam_data[,-c(1,2)]



####Load coordinates and calculate a distance matrix
coordinate <- read.csv("Documents/GitHub/Dust_project/data/metadata/points_coordinate.csv")
coordinate=coordinate[order(coordinate$id),]
rownames(coordinate)=coordinate$id
coordinate=coordinate[-26,]
coordinates(coordinate) <- cbind(coordinate$X , coordinate$Y) #### Transform it into a spatial file
proj4string(coordinate) = CRS("+proj=longlat +datum=WGS84 +no_defs")
distmatrix<-geosphere::distm(coordinate) 

beta_div = vegdist(t(asvtab_m), method="bray")
micro.dist = as.matrix(vegdist(t(asvtab_m), method="bray"))
databac <- data.frame(Dissimilarity = as.vector(micro.dist),
                      Distance = as.vector(distmatrix))

databac = databac[which (databac$Dissimilarity!=0),]
cor.test(databac$Dissimilarity, databac$Distance)



#mantel R
beta_div = vegdist(t(asvtab_m), method="bray")
dist.mat <- as.dist(distmatrix)

mantel.rtest(beta_div, dist.mat, nrepet = 9999)




a = ggplot(databac, aes(Distance,Dissimilarity))+
  geom_point(shape = 1) +
  annotate("text", x=140000, y=0.65, label= "Mantel r = 0.12; P = 0.03") +
  stat_smooth(method = "lm", formula = y ~ x,   se = FALSE)+
  theme_bw() +
  ggtitle("Fungi") +
  #stat_cor(method = "pearson", label.x = 105000, label.y = 0.49) +
  theme_classic()
a

ggsave( "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/dist_dec_fun.pdf",
        device = "pdf", 
        a)



saveRDS(a, "/Users/gabri/OneDrive - University of Arizona/dust_project/IMAGE_ELABORATION/cor_fun.RDS")





