library(metagMisc)
library(adespatial)
library(sp)
library(HH)
library(vegan)







####increasing in dust vs not decreasing in dust
geodata  = readRDS("Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS")










#### Load the data for soil###
ASVtab_soil=readRDS("Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soilITS_ASVnames.RDS")
sam_data_soil = data.frame(ASVtab_soil@sam_data) #extract soil data from phyloseq object
ASVtab_soil_merged <- merge_samples(ASVtab_soil, sample_data(ASVtab_soil)$point, fun=mean) ####merge togheter samples at each location
ASVtab_soil_merged = phyloseq_transform_css(ASVtab_soil_merged)
asvtab_soil_merged = data.frame(t(ASVtab_soil_merged@otu_table))  #extract ASV table merged










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






high_abundance = asvtab_soil_merged[,colnames(asvtab_soil_merged) %in% rownames(high_abundance_ASVs)]
high_occupancy = asvtab_soil_merged[,colnames(asvtab_soil_merged) %in% rownames(high_occupancy_ASVs)]






####first run the models with the the whole ASV table
braycurtis=vegdist(asvtab_soil_merged, method = "bray") #calculate the distance
sam_data_m=data.frame(ASVtab_soil_merged@sam_data)
sam_data_m = sam_data_m[,-c(1,2,5,6,7,8)]
#################redo VIF since I merged the variables###
sort(vif(sam_data_m))
sort(vif(sam_data_m[,-4]), decreasing = T) # eliminate total C
sort(vif(sam_data_m[,-c(3,4)]), decreasing = T) # eliminate total N (correlates  to organic carbon)
###ok, same as before
sam_data_m = sam_data_m[,-c(3,4)]
#standardized_table
sam_data_m<-decostand(sam_data_m, method="standardize", na.rm=TRUE)
dbrda_0<-dbrda( braycurtis ~ 1, sam_data_m)  # Model with intercept only
dbrda_1<-dbrda( braycurtis ~ ., sam_data_m)  # Model with all explanatory variables
ordistep(dbrda_0, scope=formula(dbrda_1), perm.max=999, direction="forward") #forward selection
db_RDA_bac_env=dbrda( braycurtis ~ pH + Conductivity + Inorganic_carbon  , sam_data_m) ###formula of significant environmental variables
db_RDA_bac_env
set.seed(7)
anova1<-anova(db_RDA_bac_env, by="margin", permutation = 9999)
anova1
RsquareAdj(db_RDA_bac_env)
env_sig = sam_data_m[, c(1,5,4)]




#####Space#####
coordinate <- read.csv("Documents/GitHub/Dust_project/data/metadata/points_coordinate.csv")
coordinate=coordinate[order(coordinate$id),]
coordinate=coordinate[-26,]
rownames(coordinate)=coordinate$X.SampleID
coordinates(coordinate) <- cbind(coordinate$X , coordinate$Y) #### Transform it into a spatial file
proj4string(coordinate) = CRS("+proj=longlat +datum=WGS84 +no_defs")
coordinates=data.frame(coordinate)[,c(1,2)]
mem.db=dbmem(coordinates, silent = FALSE)
detrended = dbrda(braycurtis~ X + Y, data = coordinate)
anova(detrended, by = "margin")
detrended  #significant trend 
detrended = resid(lm(as.matrix(braycurtis)~ X + Y, data = coordinate))

#Model     2   0.9619 1.4297  0.004 **
#detrended = resid(lm(as.matrix(braycurtis)~ X + Y, data = coordinate))
set.seed(2000)
dbrda_0<-dbrda( detrended ~ 1, mem.db)  # Model with intercept only
dbrda_1<-dbrda( detrended ~ ., mem.db)  # Model with all explanatory variables
ordistep(dbrda_0, scope=formula(dbrda_1), perm.max=999, direction="forward")
sig_mem = mem.db[,c(7,12)]
#mem.db$MEM9
dbrdabac_space=dbrda( braycurtis ~ MEM7  , mem.db) ###formula of significant environmental variables
anova(dbrdabac_space, by = "margin") # it explains little variation but has a significant effect
RsquareAdj(dbrdabac_space)#2 percent
vp1<-varpart(braycurtis,env_sig ,sig_mem, coordinates)
plot(vp1) #### 7% env, 1 % shared, 6 % space












####first do high abundance########
high_abundance  = high_abundance[rowSums(high_abundance[])>0,] # eliminate rows with all zeros
braycurtis=vegdist(high_abundance, method = "bray") #calculate the distance
sam_data_m=data.frame(ASVtab_soil_merged@sam_data)
sam_data_m = sam_data_m[,-c(1,2,5,6,7,8)]
###additional step for 1 threshold, eliminate rows with that are not in the other file
sam_data_m = sam_data_m[-c(12),]


#################redo VIF since I merged the variables###
sort(vif(sam_data_m))
sort(vif(sam_data_m[,-4]), decreasing = T) # eliminate total C
sort(vif(sam_data_m[,-c(3,4)]), decreasing = T) # eliminate total N (correlates  to organic carbon)
###ok, same as before
sam_data_m = sam_data_m[,-c(3,4)]
#standardized_table
sam_data_m<-decostand(sam_data_m, method="standardize", na.rm=TRUE)
dbrda_0<-dbrda( braycurtis ~ 1, sam_data_m)  # Model with intercept only
dbrda_1<-dbrda( braycurtis ~ ., sam_data_m)  # Model with all explanatory variables
ordistep(dbrda_0, scope=formula(dbrda_1), perm.max=999, direction="forward") #forward selection
db_RDA_bac_env=dbrda( braycurtis ~ pH  , sam_data_m) ###formula of significant environmental variables
db_RDA_bac_env
set.seed(7)
anova1<-anova(db_RDA_bac_env, by="margin", permutation = 9999)
anova1
RsquareAdj(db_RDA_bac_env)
env_sig = sam_data_m[, c(1)]






#####Space#####
coordinate <- read.csv("Documents/GitHub/Dust_project/data/metadata/points_coordinate.csv")
coordinate=coordinate[order(coordinate$id),]
coordinate=coordinate[-c(26,12),]
#coordinate= coordinate[-c(12),] #additional step for 0.5 threshold
rownames(coordinate)=coordinate$X.SampleID
coordinates(coordinate) <- cbind(coordinate$X , coordinate$Y) #### Transform it into a spatial file
proj4string(coordinate) = CRS("+proj=longlat +datum=WGS84 +no_defs")
coordinates=data.frame(coordinate)[,c(1,2)]
mem.db=dbmem(coordinates, silent = FALSE)
set.seed(123)
detrended = dbrda(braycurtis~ X + Y, data = coordinate)
anova(detrended, by = "margin")
detrended #no significant trend 
#detrended = resid(lm(as.matrix(braycurtis)~ X + Y, data = coordinate))

set.seed(2000)
dbrda_0<-dbrda( braycurtis ~ 1, mem.db)  # Model with intercept only
dbrda_1<-dbrda( braycurtis ~ ., mem.db)  # Model with all explanatory variables
ordistep(dbrda_0, scope=formula(dbrda_1), perm.max=999, direction="forward")
sig_mem = mem.db[,c(3,11)]
#mem.db$MEM9
#nothing significnat
dbrdabac_space=dbrda( braycurtis ~ MEM3 + MEM11, mem.db) ###formula of significant environmental variables
anova(dbrdabac_space, by = "margin") # not significant
RsquareAdj(dbrdabac_space)#4 percent
vp1<-varpart(braycurtis,env_sig ,sig_mem)
vp1<-varpart(braycurtis,env_sig ,sig_mem)
plot(vp1) 

# environment 4 % / space  5 %


















##### High occupancy
braycurtis=vegdist(high_occupancy, method = "bray") #calculate the distance
sam_data_m=data.frame(ASVtab_soil_merged@sam_data)
sam_data_m = sam_data_m[,-c(1,2,5,6,7,8)]
#################redo VIF since I merged the variables###
sort(vif(sam_data_m))
sort(vif(soil[,-4]), decreasing = T) # eliminate total C
sort(vif(soil[,-c(3,4)]), decreasing = T) # eliminate total N (correlates  to organic carbon)
###ok, same as before
sam_data_m = sam_data_m[,-c(3,4)]
#standardized_table
sam_data_m<-decostand(sam_data_m, method="standardize", na.rm=TRUE)
dbrda_0<-dbrda( braycurtis ~ 1, sam_data_m)  # Model with intercept only
dbrda_1<-dbrda( braycurtis ~ ., sam_data_m)  # Model with all explanatory variables
ordistep(dbrda_0, scope=formula(dbrda_1), perm.max=999, direction="forward") #forward selection
db_RDA_bac_env=dbrda( braycurtis ~ pH + C.N  , sam_data_m) ###formula of significant environmental variables
db_RDA_bac_env
set.seed(7)
anova1<-anova(db_RDA_bac_env, by="margin", permutation = 9999)
anova1
RsquareAdj(db_RDA_bac_env)
env_sig = sam_data_m[, c(1, 3)]




library(adespatial)
library(sp)
#####Space#####
coordinate <- read.csv("Documents/GitHub/Dust_project/data/metadata/points_coordinate.csv")
coordinate=coordinate[order(coordinate$id),]
coordinate=coordinate[-26,]
rownames(coordinate)=coordinate$X.SampleID
coordinates(coordinate) <- cbind(coordinate$X , coordinate$Y) #### Transform it into a spatial file
proj4string(coordinate) = CRS("+proj=longlat +datum=WGS84 +no_defs")
coordinates=data.frame(coordinate)[,c(1,2)]
mem.db=dbmem(coordinates, silent = FALSE)
detrended = dbrda(braycurtis~ X + Y, data = coordinate)
anova(detrended, by = "margin")
detrended #significant trend 
detrended = resid(lm(as.matrix(braycurtis)~ X + Y, data = coordinate))
#Model     2   0.9619 1.4297  0.004 **
#detrended = resid(lm(as.matrix(braycurtis)~ X + Y, data = coordinate))
set.seed(2000)
dbrda_0<-dbrda( detrended ~ 1, mem.db)  # Model with intercept only
dbrda_1<-dbrda( detrended ~ ., mem.db)  # Model with all explanatory variables
ordistep(dbrda_0, scope=formula(dbrda_1), perm.max=999, direction="forward")
sig_mem = mem.db[,c(9,12,10)]
#mem.db$MEM9
set.seed(125)
dbrdabac_space=dbrda( detrended ~  MEM9 +MEM12 + MEM10 , mem.db) ###formula of significant environmental variables
anova(dbrdabac_space, by = "margin") # it explains little variation but has a significant effect
RsquareAdj(dbrdabac_space)#2 percent
vp1<-varpart(braycurtis,env_sig ,sig_mem, coordinates)
plot(vp1) 


### 11% environment, 3 % shared, 18 % space



###########barplots
barplots = data.frame(matrix(ncol = 3, nrow = 3))
colnames(barplots) = c("Soil parameters","Soil/spatial","Spatial variables")
rownames(barplots) = c("All ASVs", "Generalists", "Specialists" )
#insert values manually ( I know, lame)
matrix = rbind(c(7,1,6 ),c(11,3,18),c(4,0,5))
barplots[1:3,1:3] = matrix

library(reshape2)
barplots$id = rownames(barplots)
melted = melt(barplots)

library(RColorBrewer)
library(ggplot2)
melted$id = factor(melted$id, levels = c("All ASVs", "Generalists", "Specialists") )



space_analysis_fungi = ggplot(melted, aes(fill=variable, y=value, x=id)) + 
  geom_bar(position="stack", stat="identity")  +theme_classic() + ylab("Explained variance (%)") + xlab(NULL) +
  theme(legend.position="bottom", legend.title = element_blank()) +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_brewer(palette = "Dark2") +ggtitle("Fungi")
space_analysis_fungi

ggsave(filename = "Expalined_variance_fungi.PDF", path = "/Users/gabri/OneDrive - University of Arizona/dust_project/Second_draft/new_figures/" , device = "pdf")



##change into spatial variables
#













