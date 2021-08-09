##########Script geographical values
library(vegan)
library(sp)
library(adehabitatHR) 
library(scales) # Helps make polygons partly transparent using the alpha argument below
library(geosphere)
library(tidyverse)





########ABUNDANCE#################################
#### calculate average Abundance per point
asvtab<-readRDS("Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soil16sASVnames.RDS")
asvtable_samples = data.frame(t(asvtab@otu_table)) ####take ASV table out of phyloseq object
#### create a dataframe where to store all the variables###
relative_abundances = data.frame(matrix(ncol = 3, nrow = nrow(asvtable_samples))) #create dataframe
colnames(relative_abundances) = c("Ab_bysamples","Ab_bypoint","Ab_total") #give name to the columns
rownames(relative_abundances) = rownames(asvtable_samples) # give names to the rows
#per samples
relative_abundances_samples = decostand(asvtable_samples, method = "total" , MARGIN = 2) #calculate relative abundances per sample 
relative_abundances_samples[relative_abundances_samples == 0] <- NA #create a vector to store the variables
relative_abundances$Ab_bysamples = rowMeans(relative_abundances_samples, na.rm = TRUE)  #calculate means and store it in a vector
#per point
asvtabmerged=merge_samples(asvtab, as.factor(sample_data(asvtab)$point), fun=mean) #merge samples into points (3 repetitions merged into one)
asvtabmerged <- prune_taxa(taxa_sums(asvtabmerged) > 0, asvtabmerged) #merge samples into points (3 repetitions merged into one)
asvtable_points = data.frame(t(asvtabmerged@otu_table)) #extract dataframe from  phyloseq
relative_abundances_points = decostand(asvtable_points, method = "total" , MARGIN = 2) #calculate relative abundance
relative_abundances_points[relative_abundances_points == 0] <- NA  #create a vector to store the variables
relative_abundances$Ab_bypoint = rowMeans(relative_abundances_points, na.rm = TRUE) #calculate means and store it in a vector
# total 
total_sums=sum(asvtable_samples) #sum everything
relative_abundances$Ab_total=rowSums(asvtable_samples)/total_sums #calculate the relative abundance 








#############OCCUPANCY###############################################################
#####Now make a dataframe to use for storing the data.
# I will put the occupancy for both each sample (each 3 repetition) and each point (repetitions merged)
occupancy = data.frame(matrix(ncol = 4, nrow = nrow(asvtable_samples)))
colnames(occupancy) = c("OC_samples","OC_samples_percent", "OC_points", "OC_points_percent")
rownames(occupancy) = rownames(asvtable_samples)
# Now I will calculate the occupancy per 87 samples 
occupancy$OC_samples=rowSums(asvtable_samples != 0)
occupancy$OC_samples_percent=(occupancy$OC_samples/86)*100


####now let us calculate the occupancy for merged points
occupancy$OC_points=rowSums(asvtable_points != 0)
occupancy$OC_points_percent=(occupancy$OC_points/29)*100














#############GEOGRAPHIC RANGE########################################################
###minimum convex polygon
#read the coordinate file
coordinates = read.csv("Documents/GitHub/Dust_project/data/metadata/points_coordinate.csv")
coordinates = coordinates[order(coordinates$id),]#order by ID
coordinates = coordinates[-26,]#eliminate point 26 playa
rownames(coordinates) = colnames(asvtable_points)


areas = data.frame(matrix(nrow = 0,ncol = 3))########create a dataframe where to store the variables
colnames(areas)=c("ID","X","Y")########name the columns
asvtab5= asvtable_points[which(rowSums(asvtable_points != 0) > 5) , ] #####select only points with more than 5 occupancy (function requires this)


###now sum the coordinates per each point
for (q in 1:length(rownames(asvtab5)))
{
  ASV = rownames(asvtab5)[q] 
  vector_of_coordinates = t(asvtab5[ASV,])
  coordinates_point = data.frame(matrix(nrow = length(vector_of_coordinates),ncol = 3))
  colnames(coordinates_point)=c("ID","X","Y")
  for (i in 1:length(vector_of_coordinates))
  {
    if(vector_of_coordinates[i] != 0 )
    { 
      coordinates_point[i,1] = ASV
      coordinates_point[i,2] = coordinates$X[i] 
      coordinates_point[i,3] = coordinates$Y[i] 
    }
  } 
  coordinates_point <-na.omit(coordinates_point)  
  areas = rbind(areas, coordinates_point)
}




coordinates(areas) <- cbind(areas$X , areas$Y) #### Transform it into a spatial file
areas = areas[,-c(2,3)]
proj4string(areas) = CRS("+proj=longlat +datum=WGS84 +no_defs")
areas <- spTransform(areas, CRS("+init=epsg:3857"))
plot(areas)

turtles.mcp <- mcp(areas, percent = 100, unout = "km2")
turtles.mcp
plot(turtles.mcp, col = alpha(1:5, 0.5), add = FALSE)
areas1 = data.frame(turtles.mcp@data)


#######Calculate Maximum distances#######
#####select only points with more than 2 occupancy (distances requires this)
asvtab2= asvtable_points[which(rowSums(asvtable_points != 0) > 2) , ]
distances = data.frame(matrix(nrow = nrow(asvtab2),ncol = 1))
colnames(distances)=c("max_dist")
rownames(distances)=rownames(asvtab2)

for (q in 1:length(rownames(asvtab2)))
{
  ASV = rownames(asvtab2)[q] 
  vector_of_coordinates = t(asvtab2[ASV,])
  coordinates_point = data.frame(matrix(nrow = length(vector_of_coordinates),ncol = 2))
  colnames(coordinates_point)=c("X","Y")
  for (i in 1:length(vector_of_coordinates))
  {
    if(vector_of_coordinates[i] != 0 )
    { 
      coordinates_point[i,1] = coordinates$X[i] 
      coordinates_point[i,2] = coordinates$Y[i] 
    }
  } 
  coordinates_point <-na.omit(coordinates_point)  
  distances[ASV,] = max(distm(coordinates_point))
}

distances = distances/1000















########save final table#############
final = cbind(occupancy, relative_abundances)
final$distances = NA
for (i in 1:nrow(final))  {final$distances[i] = distances[rownames(final)[i], 1] }
final$geo_range = NA
for (i in 1:nrow(final))  {final$geo_range[i] = areas1[rownames(final)[i], 2] }

saveRDS(final, file = "Documents/GitHub/Dust_project/data/metadata/geo_dist_data.RDS")



#############################################################################
###########################FUNGI#############################################
##############################################################################

##########Script geographical values
library(vegan)
library(sp)
library(adehabitatHR) 
library(scales) # Helps make polygons partly transparent using the alpha argument below
library(geosphere)
library(tidyverse)





########ABUNDANCE#################################
#### calculate average Abundance per point
asvtab<-readRDS("Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soilITS_ASVnames.RDS")
asvtable_samples = data.frame(t(asvtab@otu_table)) ####take ASV table out of phyloseq object

#### create a dataframe where to store all the variables###
relative_abundances = data.frame(matrix(ncol = 3, nrow = nrow(asvtable_samples))) #create dataframe
colnames(relative_abundances) = c("Ab_bysamples","Ab_bypoint","Ab_total") #give name to the columns
rownames(relative_abundances) = rownames(asvtable_samples) # give names to the rows
#per samples
relative_abundances_samples = decostand(asvtable_samples, method = "total" , MARGIN = 2) #calculate relative abundances per sample 
relative_abundances_samples[relative_abundances_samples == 0] <- NA #create a vector to store the variables
relative_abundances$Ab_bysamples = rowMeans(relative_abundances_samples, na.rm = TRUE)  #calculate means and store it in a vector
#per point
asvtabmerged=merge_samples(asvtab, as.factor(sample_data(asvtab)$point), fun=mean) #merge samples into points (3 repetitions merged into one)
asvtabmerged <- prune_taxa(taxa_sums(asvtabmerged) > 0, asvtabmerged) #merge samples into points (3 repetitions merged into one)
asvtable_points = data.frame(t(asvtabmerged@otu_table)) #extract dataframe from  phyloseq
relative_abundances_points = decostand(asvtable_points, method = "total" , MARGIN = 2) #calculate relative abundance
relative_abundances_points[relative_abundances_points == 0] <- NA  #create a vector to store the variables
relative_abundances$Ab_bypoint = rowMeans(relative_abundances_points, na.rm = TRUE) #calculate means and store it in a vector
# total 
total_sums=sum(asvtable_samples) #sum everything
relative_abundances$Ab_total=rowSums(asvtable_samples)/total_sums #calculate the relative abundance 








#############OCCUPANCY###############################################################
#####Now make a dataframe to use for storing the data.
# I will put the occupancy for both each sample (each 3 repetition) and each point (repetitions merged)
occupancy = data.frame(matrix(ncol = 4, nrow = nrow(asvtable_samples)))
colnames(occupancy) = c("OC_samples","OC_samples_percent", "OC_points", "OC_points_percent")
rownames(occupancy) = rownames(asvtable_samples)
# Now I will calculate the occupancy per 87 samples 
occupancy$OC_samples=rowSums(asvtable_samples != 0)
occupancy$OC_samples_percent=(occupancy$OC_samples/85)*100


####now let us calculate the occupancy for merged points
occupancy$OC_points=rowSums(asvtable_points != 0)
occupancy$OC_points_percent=(occupancy$OC_points/29)*100














#############GEOGRAPHIC RANGE########################################################
###minimum convex polygon
#read the coordinate file
coordinates = read.csv("Documents/GitHub/Dust_project/data/metadata/points_coordinate.csv")
coordinates = coordinates[order(coordinates$id),]#order by ID
coordinates = coordinates[-26,]#eliminate point 26 playa
rownames(coordinates) = colnames(asvtable_points)


areas = data.frame(matrix(nrow = 0,ncol = 3))########create a dataframe where to store the variables
colnames(areas)=c("ID","X","Y")########name the columns
asvtab5= asvtable_points[which(rowSums(asvtable_points != 0) > 5) , ] #####select only points with more than 5 occupancy (function requires this)


###now sum the coordinates per each point
for (q in 1:length(rownames(asvtab5)))
{
  ASV = rownames(asvtab5)[q] 
  vector_of_coordinates = t(asvtab5[ASV,])
  coordinates_point = data.frame(matrix(nrow = length(vector_of_coordinates),ncol = 3))
  colnames(coordinates_point)=c("ID","X","Y")
  for (i in 1:length(vector_of_coordinates))
  {
    if(vector_of_coordinates[i] != 0 )
    { 
      coordinates_point[i,1] = ASV
      coordinates_point[i,2] = coordinates$X[i] 
      coordinates_point[i,3] = coordinates$Y[i] 
    }
  } 
  coordinates_point <-na.omit(coordinates_point)  
  areas = rbind(areas, coordinates_point)
}




coordinates(areas) <- cbind(areas$X , areas$Y) #### Transform it into a spatial file
areas = areas[,-c(2,3)]
proj4string(areas) = CRS("+proj=longlat +datum=WGS84 +no_defs")
areas <- spTransform(areas, CRS("+init=epsg:3857"))
plot(areas)

turtles.mcp <- mcp(areas, percent = 100, unout = "km2")
turtles.mcp
plot(turtles.mcp, col = alpha(1:5, 0.5), add = FALSE)
areas1 = data.frame(turtles.mcp@data)



#######Calculate Maximum distances#######
#####select only points with more than 2 occupancy (distances requires this)
asvtab2= asvtable_points[which(rowSums(asvtable_points != 0) > 2) , ]
distances = data.frame(matrix(nrow = nrow(asvtab2),ncol = 1))
colnames(distances)=c("max_dist")
rownames(distances)=rownames(asvtab2)

for (q in 1:length(rownames(asvtab2)))
{
  ASV = rownames(asvtab2)[q] 
  vector_of_coordinates = t(asvtab2[ASV,])
  coordinates_point = data.frame(matrix(nrow = length(vector_of_coordinates),ncol = 2))
  colnames(coordinates_point)=c("X","Y")
  for (i in 1:length(vector_of_coordinates))
  {
    if(vector_of_coordinates[i] != 0 )
    { 
      coordinates_point[i,1] = coordinates$X[i] 
      coordinates_point[i,2] = coordinates$Y[i] 
    }
  } 
  coordinates_point <-na.omit(coordinates_point)  
  distances[ASV,] = max(distm(coordinates_point))
}

distances = distances/1000















########save final table#############
final = cbind(occupancy, relative_abundances)
final$distances = NA
for (i in 1:nrow(final))  {final$distances[i] = distances[rownames(final)[i], 1] }
final$geo_range = NA
for (i in 1:nrow(final))  {final$geo_range[i] = areas1[rownames(final)[i], 2] }

saveRDS(final, file = "Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS")






