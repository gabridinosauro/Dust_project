#####Do a baplot of compisition barplots for generalists and specialists


#####test different thresholds for community composition#######
####High abundance and High occupancy idea######
geodata  = readRDS("Documents/GitHub/Dust_project/data/metadata/geo_dist_data.RDS") # to produce this file run the code "geographic_calculations.RDS"
ASVtab_soil=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soil16sASVnames.RDS")
taxatab = data.frame(ASVtab_soil@tax_table)
taxatab[is.na(taxatab)] <- 'unknown'
tax_table(ASVtab_soil) = as.matrix(taxatab)



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








###select the ASVs based on high occupancy
selection = round((nrow(geodata)/100)*threshold) ###define the number of ASVs to pick
high_occupancy_ASVs = geodata[order(-geodata$OC_points_percent),] ### order the table
high_occupancy_ASVs = high_occupancy_ASVs[1:selection,] ### extract those with the highest occupancy
min = min(high_occupancy_ASVs$OC_points_percent) ####correct for those having the same occupancy
high_occupancy_ASVs = geodata[which(geodata$OC_points_percent >= min),] ##finally select them



####First make a boxplot for all the ASVs
dat.aglo = tax_glom(ASVtab_soil, taxrank = "Phylum")
dat.dataframe = psmelt(dat.aglo)
dat.dataframe$Phylum= as.character(dat.dataframe$Phylum)
###keep 10 most abundant phyla
ten<-aggregate(Abundance~Phylum,data= dat.dataframe, FUN=sum)
ten<-ten[order(ten$Abundance, decreasing=TRUE),]
ten<-as.vector(ten$Phylum[1:10])
`%nin%` = Negate(`%in%`) # create a negation, for not included (there is not in R)
for ( i  in 1:length(dat.dataframe$Phylum)) #substitute those not in the top 10
  if(dat.dataframe$Phylum[i] %nin% ten) dat.dataframe$Phylum[i]<-"others"

dat.dataframe<-aggregate(Abundance~Phylum,data= dat.dataframe, FUN=sum)
final = dat.dataframe
final$ASVs = "all"











####Now plot high occupancy
high_occupancy_ASVtab =  subset_taxa(ASVtab_soil, taxa_names(ASVtab_soil) %in% rownames(high_occupancy_ASVs))
dat.aglo = tax_glom(high_occupancy_ASVtab, taxrank = "Phylum")
dat.dataframe = psmelt(dat.aglo)
dat.dataframe$Phylum= as.character(dat.dataframe$Phylum)
for ( i  in 1:length(dat.dataframe$Phylum)) #substitute those not in the top 10
  if(dat.dataframe$Phylum[i] %nin% ten) dat.dataframe$Phylum[i]<-"others"

dat.dataframe<-aggregate(Abundance~Phylum,data= dat.dataframe, FUN=sum)
final_occ = dat.dataframe
final_occ$ASVs = "Generalists"


####Now plot high abundance
high_abundance_ASVtab =  subset_taxa(ASVtab_soil, taxa_names(ASVtab_soil) %in% rownames(high_abundance_ASVs))
dat.aglo = tax_glom(high_abundance_ASVtab, taxrank = "Phylum")
dat.dataframe = psmelt(dat.aglo)
dat.dataframe$Phylum= as.character(dat.dataframe$Phylum)
for ( i  in 1:length(dat.dataframe$Phylum)) #substitute those not in the top 10
  if(dat.dataframe$Phylum[i] %nin% ten) dat.dataframe$Phylum[i]<-"others"

dat.dataframe<-aggregate(Abundance~Phylum,data= dat.dataframe, FUN=sum)
final_ab = dat.dataframe
final_ab$ASVs = "Specialists"


final = rbind(final, final_occ, final_ab )


a<-ggplot(final, aes( x = ASVs, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity", position = "fill") +
  theme_bw() +
  scale_fill_brewer(palette="Paired")  +
  theme(axis.text.x = element_text(angle =45, hjust = 1))  + ylab("Relative abundance")
a




















#####FUNGI##############################################################################################
####High abundance and High occupancy idea######
geodata  = readRDS("Documents/GitHub/Dust_project/data/metadata/geo_dist_dataITS.RDS") # to produce this file run the code "geographic_calculations.RDS"
ASVtab_soil=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soilITS_ASVnames.RDS")
taxatab = data.frame(ASVtab_soil@tax_table)
taxatab[is.na(taxatab)] <- 'unknown'
tax_table(ASVtab_soil) = as.matrix(taxatab)







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








###select the ASVs based on high occupancy
selection = round((nrow(geodata)/100)*threshold) ###define the number of ASVs to pick
high_occupancy_ASVs = geodata[order(-geodata$OC_points_percent),] ### order the table
high_occupancy_ASVs = high_occupancy_ASVs[1:selection,] ### extract those with the highest occupancy
min = min(high_occupancy_ASVs$OC_points_percent) ####correct for those having the same occupancy
high_occupancy_ASVs = geodata[which(geodata$OC_points_percent >= min),] ##finally select them



####First make a boxplot for all the ASVs
dat.aglo = tax_glom(ASVtab_soil, taxrank = "Class")
dat.dataframe = psmelt(dat.aglo)
dat.dataframe$Class= as.character(dat.dataframe$Class)
###keep 10 most abundant phyla
ten<-aggregate(Abundance~Class,data= dat.dataframe, FUN=sum)
ten<-ten[order(ten$Abundance, decreasing=TRUE),]
ten<-as.vector(ten$Class[1:10])
`%nin%` = Negate(`%in%`) # create a negation, for not included (there is not in R)
for ( i  in 1:length(dat.dataframe$Class)) #substitute those not in the top 10
  if(dat.dataframe$Class[i] %nin% ten) dat.dataframe$Class[i]<-"others"

dat.dataframe<-aggregate(Abundance~Class,data= dat.dataframe, FUN=sum)
final = dat.dataframe
final$ASVs = "all"











####Now plot high occupancy
high_occupancy_ASVtab =  subset_taxa(ASVtab_soil, taxa_names(ASVtab_soil) %in% rownames(high_occupancy_ASVs))
dat.aglo = tax_glom(high_occupancy_ASVtab, taxrank = "Class")
dat.dataframe = psmelt(dat.aglo)
dat.dataframe$Class= as.character(dat.dataframe$Class)
for ( i  in 1:length(dat.dataframe$Class)) #substitute those not in the top 10
  if(dat.dataframe$Class[i] %nin% ten) dat.dataframe$Class[i]<-"others"

dat.dataframe<-aggregate(Abundance~Class,data= dat.dataframe, FUN=sum)
final_occ = dat.dataframe
final_occ$ASVs = "Generalists"


####Now plot high abundance
high_abundance_ASVtab =  subset_taxa(ASVtab_soil, taxa_names(ASVtab_soil) %in% rownames(high_abundance_ASVs))
dat.aglo = tax_glom(high_abundance_ASVtab, taxrank = "Class")
dat.dataframe = psmelt(dat.aglo)
dat.dataframe$Class= as.character(dat.dataframe$Class)
for ( i  in 1:length(dat.dataframe$Class)) #substitute those not in the top 10
  if(dat.dataframe$Class[i] %nin% ten) dat.dataframe$Class[i]<-"others"

dat.dataframe<-aggregate(Abundance~Class,data= dat.dataframe, FUN=sum)
final_ab = dat.dataframe
final_ab$ASVs = "Specialists"


final = rbind(final, final_occ, final_ab )


b<-ggplot(final, aes( x = ASVs, y=Abundance, fill=Class)) + 
  geom_bar(stat="identity", position = "fill") +
  theme_bw() +
  scale_fill_brewer(palette="Paired")  +
  theme(axis.text.x = element_text(angle =45, hjust = 1))  + ylab("Relative abundance")
b


ggarrange(a,b, labels = "AUTO")
