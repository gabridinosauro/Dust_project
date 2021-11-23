####Correlation matrix soil parameters and 
####NMDS (figure 2)
library(vegan)
library(ggplot2)
library(microbiome)
library(phyloseq)
library(metagenomeSeq)
library(HH)
library(stringr)
library(GGally)
library(ggpubr)
library(ggrepel)
library(corrplot)
library(sp)
library(adespatial)





#### Load the data for soil###
ASVtab_soil=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soil16sASVnames.RDS")
asvtab_css = metagMisc::phyloseq_transform_css(ASVtab_soil) #css transformation
asvtab_css = data.frame(t(asvtab_css@otu_table)) #extract ASV table css transformed form phyloseq object
asvtab_soil = data.frame(t(ASVtab_soil@otu_table)) #extract ASVt table not transformed
ASVtab_soil_merged <- merge_samples(ASVtab_soil, sample_data(ASVtab_soil)$point, fun=mean) ####merge togheter samples at each location
ASVtab_soil_merged = metagMisc::phyloseq_transform_css(ASVtab_soil_merged)
asvtab_soil_merged = data.frame(t(ASVtab_soil_merged@otu_table))  #extract ASV table merged
sam_data_soil = data.frame(ASVtab_soil@sam_data)











#####Correlation matrix, soil variables#####
soil = sam_data_soil[,-c(1,2,5,6,7,8)]###drop some variables I do not need (NDVI, sand, micro_tot) (variables are incorrect)
colnames(soil) = c("pH","Macro/Micro aggregates","total N","total C", "C/N ratio", "Conductivity", "Inorganic C", "Organic C") # insert better names
cormatrix=cor(soil, method="spearman", use="pairwise.complete.obs") #spearman correlation matrix
res1 <- cor.mtest(soil, conf.level = .95, use="pairwise.complete.obs") #correlation tests
pAdj <- p.adjust(c(res1[[1]]), method = "BH") #adjust p-values
resAdj <- matrix(pAdj, ncol = dim(res1[[1]])[1]) #make a matrix of adjusted p-valiues
plot1 = corrplot(cormatrix, type = "upper", method= "number", p.mat = resAdj, sig.level = .05) # plot the correlation matrix















#####calculate community indexes bacteria
Richness.rar<-specnumber(rrarefy(t(asvtab_soil), 14500))# richness rarefied 
Shannon.rar = vegan::diversity(rrarefy(t(asvtab_soil), 14500))#diversity
Evenness.rar = Shannon.rar/log(Richness.rar)##eveness
bacteria.ecolind= cbind(Richness.rar,Shannon.rar,Evenness.rar)#bid all togheter









#####calculate community indexes for fungi
#load the data
ASVtab_soil_fungi=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soilITS_ASVnames.RDS")
asvtab_css_fungi = metagMisc::phyloseq_transform_css(ASVtab_soil_fungi) #css transformation
asvtab_css_fungi = data.frame(t(asvtab_css_fungi@otu_table)) #extract ASV table css transformed form phyloseq object
asvtab_soil_fungi = data.frame(t(ASVtab_soil_fungi@otu_table)) #extract ASVt table not transformed
ASVtab_soil_merged_fungi <- merge_samples(ASVtab_soil_fungi, sample_data(ASVtab_soil_fungi)$point, fun=mean) ####merge togheter samples at each location
ASVtab_soil_merged_fungi = metagMisc::phyloseq_transform_css(ASVtab_soil_merged)
asvtab_soil_merged_fungi = data.frame(t(ASVtab_soil_merged@otu_table))  #extract ASV table merged
#calculate community indexes
Richness.rar.fung = specnumber(rrarefy(t(asvtab_soil_fungi), 12000))
Shannon.rar.fung = vegan::diversity(rrarefy(t(asvtab_soil_fungi), 12000))
Evenness.rar.fung = Shannon.rar.fung/log(Richness.rar.fung)
fungi.ecolind= cbind(Richness.rar.fung,Richness.rar.fung,Evenness.rar.fung)












###bind bacteria and fungi together
ecoldindicators = merge(bacteria.ecolind, fungi.ecolind, all = TRUE, by = 'row.names')
rownames(ecoldindicators) = ecoldindicators[,1]
ecoldindicators = ecoldindicators[,-1]














###correlation matrix between bacteria?fungi community indexes and soil parameters
A_data2 <- data.frame(ecoldindicators )               #####first  matrix####
colnames(A_data2) = c("Bacteria-archaea richness", "Bacteria-archaea diversity (H)", "Bacteria-Archaea evenness (J)", "Fungi richness", "Fungi diversity (H)", "Fungi evenness (J)")
### have use the original metadata file because bacteria and fungi have different rows
metadata_soil=read.csv("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/metadata_soil.csv", row.names = 1)
A_data3 <- metadata_soil[-c(76,77,78),-c(1,2,5,6,7,8)]         #####second matrix####


M <- cor(A_data2, A_data3, use="pairwise.complete.obs", method="spearman") ####correlation coefficient###
# significance test (fuction definition)
cor2.mtest <- function(mat1,mat2, conf.level = 0.95){
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  n <-ncol(mat1)
  l <-ncol (mat2)
  p.mat  <- matrix(NA, n, l)
  
  for(i in 1:n){
    for(j in 1:l){
      tmp <- cor.test(mat1[,i], mat2[,j], conf.level = conf.level)
      p.mat[i,j]<- tmp$p.value
    }
  }
  return(list(p.mat,nrow=n,ncol=l))
}
res1 <- cor2.mtest(A_data2,A_data3,0.95)  ######make significance test
pAdj <- p.adjust(c(res1[[1]]), method = "BH") #adjust p-values
resAdj <- matrix(pAdj, ncol = dim(res1[[1]])[1]) #make a matrix of adjusted p-valiues
plot2 = corrplot(M ,method="number",  p.mat = t(resAdj), sig.level=0.05) ####plot everything











######plot correlation matrix between community indexes
cormatrix=cor(A_data2, method="spearman", use="pairwise.complete.obs") #spearman correlation matrix
res1 <- cor.mtest(A_data2, conf.level = .95, use="pairwise.complete.obs") #correlation tests
pAdj <- p.adjust(c(res1[[1]]), method = "BH") #adjust p-values
resAdj <- matrix(pAdj, ncol = dim(res1[[1]])[1]) #make a matrix of adjusted p-valiues
plot1 = corrplot(cormatrix, type = "upper", method= "number", p.mat = resAdj, sig.level = .05) # plot the correlation matrix

a = ggplot(A_data2, aes(x =`Bacteria-archaea richness`, y =`Fungi richness` )) +geom_point() +geom_smooth(span = 2) +theme_classic()
b = ggplot(A_data2, aes(x =`Bacteria-archaea diversity (H)`, y =`Fungi diversity (H)` )) +geom_point() +geom_smooth(span = 2) +theme_classic()
par(mfrow=c(1,1))
ggarrange(a,b, labels = "AUTO")















##################VIF SOIl
sort(vif(soil))
sort(vif(soil[,-4]), decreasing = T) # eliminate total C
sort(vif(soil[,-c(3,4)]), decreasing = T) # eliminate total N (correlates  to organic carbon)



###all values are lower than 5 (we do not want valued higher than 5)
sam_data = soil[,-c(3,4)]
sam_data_bac = decostand(sam_data, method = "standardize" )
####
bacteria_bc = vegdist(asvtab_css, method="bray") ###somehow it does not like it
permanovabac = adonis(asvtab_css ~ ., data= sam_data_bac) ###anyway in the function adonis I can use the method bray too
permanovabac
#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#                           Df  SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  pH                        1    3.5646  3.5646 14.1945 0.13659  0.001 ***
#  `Macro/Micro aggregates`  1    0.5366  0.5366  2.1367 0.02056  0.007 ** 
#  `C/N ratio`               1    0.6377  0.6377  2.5392 0.02443  0.003 ** 
#  Conductivity              1    0.5419  0.5419  2.1579 0.02077  0.004 ** 
#  `Inorganic C`             1    0.5725  0.5725  2.2796 0.02194  0.002 ** 
#  `Organic C`               1    0.4041  0.4041  1.6091 0.01548  0.035 *  
#  Residuals                79   19.8390  0.2511         0.76022           
#  Total                    85   26.0963                 1.00000           
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

















#####db_RDA, to avoid inflation given by the tree replicates, I used merged points (values are merged for each location togheter)
#forming a unique sample)
braycurtis=vegdist(asvtab_soil_merged, method = "bray") #calculate the distance
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
db_RDA_bac_env=dbrda( braycurtis ~ pH + Inorganic_carbon + Organic_carbon...  , sam_data_m) ###formula of significant environmental variables
db_RDA_bac_env
set.seed(7)
anova1<-anova(db_RDA_bac_env, by="margin", permutation = 9999)
anova1

RsquareAdj(db_RDA_bac_env)
env_sig = sam_data_m[, c(1,5,6)]
















#####space#####
coordinate <- read.csv("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/points_coordinate.csv")
coordinate=coordinate[order(coordinate$id),]
coordinate=coordinate[-26,]
rownames(coordinate)=coordinate$X.SampleID
coordinates(coordinate) <- cbind(coordinate$X , coordinate$Y) #### Transform it into a spatial file
proj4string(coordinate) = CRS("+proj=longlat +datum=WGS84 +no_defs")
coordinates=data.frame(coordinate)[,c(1,2)]
mem.db=dbmem(coordinates, silent = FALSE)
#writeRDS to use in another script
#saveRDS(mem.db, "/Users/gabri/OneDrive - University of Arizona/dust_project/data/soil/membac.RDS")
#test_trend 
detrended = anova(dbrda(braycurtis~ X + Y, data = coordinate))
detrended #no significant trend




set.seed(2000)
dbrda_0<-dbrda( braycurtis ~ 1, mem.db)  # Model with intercept only
dbrda_1<-dbrda( braycurtis ~ ., mem.db)  # Model with all explanatory variables
ordistep(dbrda_0, scope=formula(dbrda_1), perm.max=999, direction="forward") ####nothing
sig_mem = NA
#Variation_partitioning
#vp1<-varpart(braycurtis,env_sig ,sig_mem)
#plot(vp1)
















#Bacteria####
micro.dist = vegdist(asvtab_css, method="bray")

# Now calculate NMDS axes
micro.nmds = metaMDS(micro.dist, k=2)

# Add NMDS axes to metadata
md=data.frame(ASVtab_soil@sam_data)
md$Axis01 = micro.nmds$points[,1]
md$Axis02 = micro.nmds$points[,2]
md$richness_tot = Richness.rar

temp <- c()

for (i in unique(md$point)) {
  temp1 <- subset(md, md$point == i)
  temp2 <- aggregate(cbind(Axis01, Axis02, richness_tot) ~ point, data = temp1, mean)
  temp3 <- aggregate(cbind(Axis01, Axis02) ~ point, data = temp1, sd)
  names(temp3)[c(2,3)] <- c("NMDS1.sd", "NMDS2.sd")
  temp1 <- cbind(temp2, temp3)
  temp1$Site <- i
  temp <- rbind(temp, temp1)
}

graph = temp
rm(temp, temp1, temp2, temp3, i)
graph$pH = ASVtab_soil_merged@sam_data$pH

# Min and Max of NMDS1 and NMDS2
graph$NMDS1.min <- graph$Axis01 - graph$NMDS1.sd
graph$NMDS1.max <- graph$Axis01 + graph$NMDS1.sd
graph$NMDS2.min <- graph$Axis02 - graph$NMDS2.sd
graph$NMDS2.max <- graph$Axis02 + graph$NMDS2.sd


# Remove duplicate columns
graph <- graph[,-1]

# plot
mid<-mean(graph$pH)
stress= paste("stress", round(micro.nmds$stress, digits = 2))
bac.NMDS.plot <- ggplot(data = graph, aes(x = Axis01, y = Axis02, color = pH))+
  geom_point(aes(size = richness_tot))+
  geom_errorbar(aes(ymin = NMDS2.min, ymax = NMDS2.max))+
  geom_errorbarh(aes(xmin = NMDS1.min, xmax = NMDS1.max))+
  #geom_label_repel(aes(label = datenum), size = 3)+
  theme_classic(base_size = 15)+
  ggtitle("Bacteria/Archaea, soil samples")+
  annotate("text",x=max(md$Axis01),y=min(md$Axis02),hjust=1,label= stress,  size = 3)+
  scale_color_viridis_c(direction=-1)+
  labs(size="Richness")+
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
  
bac.NMDS.plot

























###############Fungi#################################
###all values are lower than 5 (we do not want valued higher than 5)
ASVtab_soil_fungi=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soilITS_ASVnames.RDS")
asvtab_css_fungi = metagMisc::phyloseq_transform_css(ASVtab_soil_fungi) #css transformation
asvtab_css_fungi = data.frame(t(asvtab_css_fungi@otu_table)) #extract ASV table css transformed form phyloseq object
asvtab_soil_fungi = data.frame(t(ASVtab_soil_fungi@otu_table)) #extract ASVt table not transformed
ASVtab_soil_fungi_merged <- merge_samples(ASVtab_soil_fungi, sample_data(ASVtab_soil_fungi)$point, fun=mean) ####merge togheter samples at each location
ASVtab_soil_fungi_merged = metagMisc::phyloseq_transform_css(ASVtab_soil_fungi_merged)
asvtab_soil_fungi_merged = data.frame(t(ASVtab_soil_fungi_merged@otu_table))  #extract ASV table merged
sam_data = data.frame(ASVtab_soil_fungi@sam_data)





sam_data = sam_data[,-c(1,2,5,6,7,8,9,10)]
sam_data = decostand(sam_data, method = "standardize" )
####
fungi_bc = vegdist(asvtab_css_fungi, method="bray") ###somehow it does not like it
permanovafun = adonis(asvtab_css_fungi ~ ., data= sam_data) ###anyway in the function adonis I can use the method bray too
permanovafun
#Permutation: free
#Number of permutations: 999

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#pH                 1    1.7704 1.77043  5.5002 0.05829  0.001 ***
#  macro_micro_ratio  1    0.6867 0.68675  2.1335 0.02261  0.001 ***
#  C.N                1    0.8027 0.80273  2.4938 0.02643  0.001 ***
#  Conductivity       1    0.6654 0.66535  2.0670 0.02191  0.001 ***
#  Inorganic_carbon   1    0.8101 0.81010  2.5167 0.02667  0.001 ***
#  Organic_carbon...  1    0.5295 0.52947  1.6449 0.01743  0.002 ** 
#  Residuals         78   25.1071 0.32189         0.82665           
# Total             84   30.3720                 1.00000           
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
















#####db_RDA, to avoid inflation given by the tree replicates, I used merged points (values are merged for each location togheter)
#forming a unique sample)
braycurtis=vegdist(asvtab_soil_fungi_merged, method = "bray") #calculate the distance
sam_data_m=data.frame(ASVtab_soil_fungi_merged@sam_data)
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
db_RDA_bac_env=dbrda( braycurtis ~ pH + Inorganic_carbon +  Conductivity , sam_data_m) ###formula of significant environmental variables
db_RDA_bac_env
set.seed(7)
anova1<-anova(db_RDA_bac_env, by="margin", permutation = 9999)
anova1
RsquareAdj(db_RDA_bac_env)
env_sig = sam_data_m[, c(1,4,5)] ###change it in the paper














#####Space#####
coordinate <- read.csv("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/points_coordinate.csv")
coordinate=coordinate[order(coordinate$id),]
coordinate=coordinate[-26,]
rownames(coordinate)=coordinate$X.SampleID
coordinates(coordinate) <- cbind(coordinate$X , coordinate$Y) #### Transform it into a spatial file
proj4string(coordinate) = CRS("+proj=longlat +datum=WGS84 +no_defs")
coordinates=data.frame(coordinate)[,c(1,2)]
mem.db=dbmem(coordinates, silent = FALSE)
detrended = anova(dbrda(braycurtis~ X + Y, data = coordinate))
detrended #significant trend 
#Model     2   0.9619 1.4297  0.004 **
detrended = resid(lm(as.matrix(braycurtis)~ X + Y, data = coordinate))
set.seed(2000)
dbrda_0<-dbrda( detrended ~ 1, mem.db)  # Model with intercept only
dbrda_1<-dbrda( detrended ~ ., mem.db)  # Model with all explanatory variables
ordistep(dbrda_0, scope=formula(dbrda_1), perm.max=999, direction="forward")
sig_mem = mem.db[,c(7, 12)]

dbrdabac_space=dbrda( braycurtis ~ MEM7, mem.db) ###formula of significant environmental variables
anova(dbrdabac_space, by = "margin") # it explains little variation but has a significant effect
RsquareAdj(dbrdabac_space)#2 percent
vp1<-varpart(braycurtis,env_sig ,sig_mem, coordinates)
plot(vp1)








#NMDS ordination fungi
# First calculate distance matrix with Bray-Curtis
# I need first to transform to relative abundances, as suggested by  https://doi.org/10.1371/journal.pcbi.1003531
micro.dist = vegdist(asvtab_css_fungi, method="bray")

# Now calculate NMDS axes
micro.nmds = metaMDS(micro.dist, k=2)

# Add NMDS axes to metadata
md=data.frame(ASVtab_soil_fungi@sam_data)
md$Axis01 = micro.nmds$points[,1]
md$Axis02 = micro.nmds$points[,2]
md$richness_tot = Richness.rar.fung

temp <- c()

for (i in unique(md$point)) {
  temp1 <- subset(md, md$point == i)
  temp2 <- aggregate(cbind(Axis01, Axis02, richness_tot) ~ point, data = temp1, mean)
  temp3 <- aggregate(cbind(Axis01, Axis02) ~ point, data = temp1, sd)
  names(temp3)[c(2,3)] <- c("NMDS1.sd", "NMDS2.sd")
  temp1 <- cbind(temp2, temp3)
  temp1$Site <- i
  temp <- rbind(temp, temp1)
}

graph = temp
rm(temp, temp1, temp2, temp3, i)
graph$pH = ASVtab_soil_merged_fungi @sam_data$pH

# Min and Max of NMDS1 and NMDS2
graph$NMDS1.min <- graph$Axis01 - graph$NMDS1.sd
graph$NMDS1.max <- graph$Axis01 + graph$NMDS1.sd
graph$NMDS2.min <- graph$Axis02 - graph$NMDS2.sd
graph$NMDS2.max <- graph$Axis02 + graph$NMDS2.sd


# Remove duplicate columns
graph <- graph[,-1]

# plot
mid<-mean(graph$pH)
stress= paste("stress", round(micro.nmds$stress, digits = 2))
fung.NMDS.plot <- ggplot(data = graph, aes(x = Axis01, y = Axis02, color = pH))+
  geom_point(aes(size = richness_tot))+
  geom_errorbar(aes(ymin = NMDS2.min, ymax = NMDS2.max))+
  geom_errorbarh(aes(xmin = NMDS1.min, xmax = NMDS1.max))+
  #geom_label_repel(aes(label = datenum), size = 3)+
  theme_classic(base_size = 15)+
  ggtitle("Fungi, soil samples")+
  annotate("text",x=max(md$Axis01),y=min(md$Axis02),hjust=1,label= stress,  size = 3)+
  scale_color_viridis_c(direction =  -1)+
  labs(size="Richness")+
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))

fung.NMDS.plot







rm(list=setdiff(ls(), c("fung.NMDS.plot", "bac.NMDS.plot")))









#Bacteria Dust####
ASVtab_dust_bac<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
sam_data = data.frame(ASVtab_dust_bac@sam_data)
sam_data$sample <- paste(sam_data$Point_type,sam_data$Repetition, sep = "_") ####add a column with the repetition

#####insert new asvtab into the phyloseq object
ASVtab_dust_bac@sam_data =sample_data(sam_data) 
ASVtab_dust_bac@sam_data$Repetition=factor(ASVtab_dust_bac@sam_data$Repetition)

asvtab_dust_bac = data.frame(ASVtab_dust_bac@otu_table)#extract ASV table
ASVtab_dust_bac_css = metagMisc::phyloseq_transform_css(ASVtab_dust_bac)#CSS transform
asvtab_dust_bac_css = data.frame(t(ASVtab_dust_bac_css@otu_table)) #extract CSS transformed ASV table
sam_data$total <- paste(sam_data$point,sam_data$Dust_Type)
ASVtab_dust_merged <- merge_samples(ASVtab_dust_bac, sam_data$total, fun=mean) ####merge togheter samples at each location
ASVtab_dust_merged = metagMisc::phyloseq_transform_css(ASVtab_dust_merged)
asvtab_dust_merged = data.frame(t(ASVtab_dust_merged@otu_table))  #extract ASV table merged




richness = specnumber(rrarefy(asvtab_dust_bac, 14500))


micro.dist = vegdist(asvtab_dust_bac_css, method="bray")
sample_metadata = data.frame(ASVtab_dust_bac@sam_data)
sample_metadata$total <- paste(sample_metadata$point,sample_metadata$Dust_Type)
sample_data(ASVtab_dust_bac)= sample_metadata


# Now calculate NMDS axes
micro.nmds = metaMDS(micro.dist, k=2)

# Add NMDS axes to metadata
md=sample_metadata
md$Axis01 = micro.nmds$points[,1]
md$Axis02 = micro.nmds$points[,2]
md$richness_tot = richness

temp <- c()
for (i in unique(md$total)) {
  temp1 <- subset(md, md$total == i)
  temp2 <- aggregate(cbind(Axis01, Axis02, richness_tot) ~ point, data = temp1, mean)
  temp3 <- aggregate(cbind(Axis01, Axis02) ~ point, data = temp1, sd)
  names(temp3)[c(2,3)] <- c("NMDS1.sd", "NMDS2.sd")
  temp1 <- cbind(temp2, temp3)
  temp1$Site <- i
  temp <- rbind(temp, temp1)
}

graph = temp
rm(temp, temp1, temp2, temp3, i)

graph$total=graph$Site
# Min and Max of NMDS1 and NMDS2
graph$NMDS1.min <- graph$Axis01 - graph$NMDS1.sd
graph$NMDS1.max <- graph$Axis01 + graph$NMDS1.sd
graph$NMDS2.min <- graph$Axis02 - graph$NMDS2.sd
graph$NMDS2.max <- graph$Axis02 + graph$NMDS2.sd
# Remove duplicate columns



# plot
graph1 = str_split_fixed(graph$Site, " ", 2)
colnames(graph1) = c("point","Sample type")
graph = cbind(graph,graph1)
graph <- graph[,-c(1,14)]
replacements <- data.table::data.table(
  old = c("1", "10", "11", "18"),
  new = c("PP", "SCR", "DS1", "DS2")
)

graph$Site <- replacements$new[match(graph$point, replacements$old)]

replacements <- data.table::data.table(
  old = c("Soil", "25", "75", "500", "Dust_Gen"),
  new = c("Soil", "Fine SS", "Medium SS", "Coarse SS", "WT")
)

graph$`Sample type` <- replacements$new[match(graph$`Sample type`, replacements$old)]


stress= paste("stress", round(micro.nmds$stress, digits = 2))
library(ggrepel)
bac.dust.NMDS.plot <- ggplot(data = graph, aes(x = Axis01, y = Axis02, color = Site))+
  geom_point(aes(size = richness_tot))+
  geom_errorbar(aes(ymin = NMDS2.min, ymax = NMDS2.max))+
  geom_errorbarh(aes(xmin = NMDS1.min, xmax = NMDS1.max))+
  #geom_label_repel(aes(label = datenum), size = 3)+
  theme_classic(base_size = 15)+
  ggtitle("Bacteria/Archaea, WT and SS samples ")+
  annotate("text",x=max(md$Axis01),y=min(md$Axis02),hjust=1,label= stress,  size = 3) + 
  geom_text_repel(aes(label =`Sample type`), size = 3.5, color = "black") +
  scale_colour_manual(values = c( "#52854C" ,"#E7B800", "#D16103","#00AFBB")) +
  labs(size="Richness") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))


bac.dust.NMDS.plot














#########FUNGI####################
dust=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
ASVtabt=data.frame(t(dust@otu_table))







#ASVtabt=t(ASVtab)
# convert OTU table into package format
metaSeqObject= newMRexperiment(ASVtabt) 
# CSS normalization
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
# convert CSS normalized data into data.frame-formatted OTU table (log transformed data)
OTU_read_count_CSS = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))
#ASVtabhel<-decostand(ASVtab, method="hellinger", na.rm=TRUE)
ASVtabCSS=t(OTU_read_count_CSS)
micro.dist = vegdist(ASVtabCSS, method="bray")
sample_metadata = data.frame(dust@sam_data)
sample_metadata$total <- paste(sample_metadata$point,sample_metadata$Dust_Type)
sample_data(dust)= sample_metadata

micro.dist = vegdist(ASVtabCSS, method="bray")

# Now calculate NMDS axes
micro.nmds = metaMDS(micro.dist, k=2)

# Add NMDS axes to metadata
md=data.frame(dust@sam_data)
md$Axis01 = micro.nmds$points[,1]
md$Axis02 = micro.nmds$points[,2]
md$richness_tot = specnumber(ASVtabCSS)

temp <- c()
for (i in unique(md$total)) {
  temp1 <- subset(md, md$total == i)
  temp2 <- aggregate(cbind(Axis01, Axis02, richness_tot) ~ point, data = temp1, mean)
  temp3 <- aggregate(cbind(Axis01, Axis02) ~ point, data = temp1, sd)
  names(temp3)[c(2,3)] <- c("NMDS1.sd", "NMDS2.sd")
  temp1 <- cbind(temp2, temp3)
  temp1$Site <- i
  temp <- rbind(temp, temp1)
}

graph = temp
rm(temp, temp1, temp2, temp3, i)
asvtabmerged=merge_samples(dust, as.factor(sample_data(dust)$total), fun=mean)
graph$total = asvtabmerged@sam_data$total
# Min and Max of NMDS1 and NMDS2
graph$NMDS1.min <- graph$Axis01 - graph$NMDS1.sd
graph$NMDS1.max <- graph$Axis01 + graph$NMDS1.sd
graph$NMDS2.min <- graph$Axis02 - graph$NMDS2.sd
graph$NMDS2.max <- graph$Axis02 + graph$NMDS2.sd
# Remove duplicate columns

# plot
library(stringr)
graph1 = str_split_fixed(graph$Site, " ", 2)
colnames(graph1) = c("point","Sample type")
graph = cbind(graph,graph1)
graph <- graph[,-c(1,14)]
replacements <- data.table::data.table(
  old = c("1", "10", "11", "18"),
  new = c("PP", "SCR", "DS1", "DS2")
)

graph$Site <- replacements$new[match(graph$point, replacements$old)]

replacements <- data.table::data.table(
  old = c("Soil", "25", "75", "500", "Dust_Gen"),
  new = c("Soil", "Fine SS", "Medium SS", "Coarse SS", "WT")
)

graph$`Sample type` <- replacements$new[match(graph$`Sample type`, replacements$old)]


stress= paste("stress", round(micro.nmds$stress, digits = 2))

fung.dust.NMDS.plot <- ggplot(data = graph, aes(x = Axis01, y = Axis02, color = Site))+
  geom_point(aes(size = richness_tot))+
  geom_errorbar(aes(ymin = NMDS2.min, ymax = NMDS2.max))+
  geom_errorbarh(aes(xmin = NMDS1.min, xmax = NMDS1.max))+
  #geom_label_repel(aes(label = datenum), size = 3)+
  theme_classic(base_size = 15)+
  ggtitle("Fungi, WT and SS samples")+
  annotate("text",x=max(md$Axis01),y=min(md$Axis02),hjust=1,label= stress,  size = 3) + 
  geom_text_repel(aes(label =`Sample type`), size = 3.5, color = "black") +
  scale_colour_manual(values = c( "#52854C" ,"#E7B800", "#D16103","#00AFBB")) +
  labs(size="Richness") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))


fung.dust.NMDS.plot


ggarrange(bac.dust.NMDS.plot, fung.dust.NMDS.plot, bac.NMDS.plot, fung.NMDS.plot, labels = "AUTO")

#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/3rd_draft/plots/new_figures/figure2.pdf", device = "pdf", height = 7, width = 11)
