library(dplyr)
library(phyloseq)

###load bacterial data####
asvtab<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/end_dada2/OTU16sinvplates.rds") #LOAD ASV tab 
sampleinfo<-read.csv("/Users/gabri/Documents/GitHub/Dust_project/data/end_dada2/predata16s.csv")#load mapping file
TAX=as.data.frame(readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/end_dada2/TAX16sinvplates.rds")) # load taxa file
sampleinfo <- sampleinfo[sampleinfo[,1] %in% as.vector(rownames(asvtab)),] ###select only points that have been analyzed
#check if they are the same
setdiff(sampleinfo$Sample.ID, as.vector(rownames(asvtab)))
setdiff(as.vector(rownames(asvtab)),sampleinfo$Sample.ID )




####Eliminate eaukaria, mitochondria and Chloroplast
fil.taxa = TAX[!TAX$Order %in% "Chloroplast",] %>% droplevels()
fil.taxa = fil.taxa[!fil.taxa$Family %in% "Mitochondria",] %>% droplevels()
fil.taxa = fil.taxa[!fil.taxa$Kingdom %in% "Eukaryota",] %>% droplevels()
fil.taxa = fil.taxa[!(is.na(fil.taxa$Kingdom)),] %>% droplevels()
# Remove from ASV table
chloro_asvs = rownames(TAX[TAX$Order %in% "Chloroplast",])
mito_asvs = rownames(TAX[TAX$Family %in% "Mitochondria",])
euky_asvs = rownames(TAX[TAX$Kingdom %in% "Eukaryota",])
na_asvs = rownames(TAX[is.na(TAX$Kingdom),])
#fill ASV
fil.asv = asvtab[,!colnames(asvtab) %in% chloro_asvs]
fil.asv = fil.asv[,!colnames(fil.asv) %in% mito_asvs]
fil.asv = fil.asv[,!colnames(fil.asv) %in% euky_asvs]
fil.asv = fil.asv[,!colnames(fil.asv) %in% na_asvs]




#now I check if the points are the same or some additional point has been cancelled
setdiff(sampleinfo$Sample.ID, rownames(fil.asv))
setdiff(rownames(fil.asv),sampleinfo$Sample.ID )




#I create the phyloseq object
identical(colnames(fil.asv),row.names(fil.taxa))
sampleinfo=sampleinfo[match(row.names(fil.asv), sampleinfo$Sample.ID),]
rownames(sampleinfo)=sampleinfo$Sample.ID
sampleinfo=sampleinfo[,-1]
OTU=otu_table(fil.asv, taxa_are_rows = FALSE)
TAX=tax_table(as.matrix(fil.taxa))
physeqmod = phyloseq(OTU, TAX, sample_data(sampleinfo))
#####save sequences and change names to ASV1 ASV2 ASV3 etc
sequences=taxa_names(physeqmod)
ASVnames=sprintf("ASV%s",seq(1:ncol(physeqmod@otu_table)))
sequencesnames=rbind(ASVnames, sequences)
#saveRDS(sequencesnames, "/Users/gabri/OneDrive - University of Arizona/dust_project/data/soil/sequences.RDS")
taxa_names(physeqmod)=ASVnames





#####remove playa point (point number 26).
#####Willcox playa is a dry lake in southern Arizona. It was sampled but given the high difference with the other sampled environments,
##### total lack of vegetation. Extremely high salinity, lack of organic carbon, very High pH
##### it will not be used since it is an outlier. 
physeqmod=prune_samples(sample_data(physeqmod)$point!=26, physeqmod)
physeqmod <- prune_taxa(taxa_sums(physeqmod) > 0, physeqmod) # Eliminate rows with only zeros





###soil phyloseq
soil<- prune_samples(sample_data(physeqmod)$Type == "Soil", physeqmod)
###insert real metadata table
metadata_soil=read.csv("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/metadata_soil.csv", row.names = 1)
metadata_soil=sample_data(metadata_soil)
soil <- merge_phyloseq(soil, metadata_soil)
soil@sam_data=soil@sam_data[,-c(2,3,4)]





###dust phyloseq
dustpoints=c(1,10,11,18)
dust<- prune_samples(sample_data(physeqmod)$point %in% dustpoints, physeqmod)
metadata_dust=read.csv("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/dust_met.csv", row.names = 1)
metadata_dust=sample_data(metadata_dust)
dust <- merge_phyloseq(dust, metadata_dust)
dust@sam_data=dust@sam_data[,-c(1,2,3)]
soil = prune_samples(sample_names(soil)!="DP24-01-S", soil) ####eliminate DP24-01-S point  - less than 10000 reads
soil <- prune_taxa(taxa_sums(soil) > 0, soil)
###summary of table
ASVtab_soil = data.frame(t(soil@otu_table))
nrow(ASVtab_soil)#16431
mean(colSums(ASVtab_soil!=0)) #number of sequences per sample 633.36
sd(colSums(ASVtab_soil!=0)) #SD 126.84
taxatab = data.frame(soil@tax_table)
ASVtab_soil$phylum = taxatab$Phylum[match(rownames(ASVtab_soil), rownames(taxatab))]
ASVtab_soil = aggregate(. ~ phylum ,data = ASVtab_soil,  FUN = sum)
rownames(ASVtab_soil) = ASVtab_soil[,1]
ASVtab_soil = ASVtab_soil[,-1]
rowsums_ASVtab_soil = rowSums(ASVtab_soil)
total_sum = sum(rowsums_ASVtab_soil)
fractions = rowsums_ASVtab_soil/total_sum
fractions = fractions[order(-fractions)]
fractions = fractions*100
fractions


###summary of table
ASVtab_dust = data.frame(t(dust@otu_table))
ASVtab_dust = ASVtab_dust[,-c(1,2,3,16,17,18,31,32,33,46,47,48)]
ASVtab_dust = ASVtab_dust[rowSums(ASVtab_dust[])>0,]
nrow(ASVtab_dust)#9448
mean(colSums(ASVtab_dust!=0)) #number of sequences per sample 786.02
sd(colSums(ASVtab_dust!=0)) #SD 176.96
taxatab = data.frame(dust@tax_table)
ASVtab_dust$phylum = taxatab$Phylum[match(rownames(ASVtab_dust), rownames(taxatab))]
ASVtab_dust = aggregate(. ~ phylum ,data = ASVtab_dust,  FUN = sum)
rownames(ASVtab_dust) = ASVtab_dust[,1]
ASVtab_dust = ASVtab_dust[,-1]
rowsums_ASVtab_dust = rowSums(ASVtab_dust)
total_sum = sum(rowsums_ASVtab_dust)
fractions = rowsums_ASVtab_dust/total_sum
fractions = fractions[order(-fractions)]
fractions = fractions*100
fractions




saveRDS(soil, "/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soil16sASVnames.RDS")
dust <- prune_taxa(taxa_sums(dust) > 0, dust)
saveRDS(dust, "/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
rm(list = ls())





















####ITS####
asvtabITS<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/end_dada2/seqtabITS.rds")
sampleinfoITS<-read.csv("/Users/gabri/Documents/GitHub/Dust_project/data/end_dada2/predata16s.csv")
TAXITS=as.data.frame(readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/end_dada2/taxaITS.rds"))
#I need to change the name of DP26 because there is a _ instead of -
rownames(asvtabITS)[rownames(asvtabITS) == "DP26_2_D500"] <- "DP26-2_D500"
rownames(asvtabITS)[rownames(asvtabITS) == "BLANKAW"] <- "BLANK"
rownames(asvtabITS)[rownames(asvtabITS) == "BLANKAW1"] <- "BLANK1"
rownames(asvtabITS)[rownames(asvtabITS) == "BLANKAW2"] <- "BLANK2"
sampleinfoITS <- sampleinfoITS[sampleinfoITS[,1] %in% as.vector(rownames(asvtabITS)),] ###select only points that have been analyzed
setdiff(sampleinfoITS$Sample.ID, as.vector(rownames(asvtabITS)))
setdiff(as.vector(rownames(asvtabITS)),sampleinfoITS$Sample.ID )
identical(colnames(asvtabITS),row.names(TAXITS))
#I create the phyloseq object
library(phyloseq)
sampleinfoITS=sampleinfoITS[match(row.names(asvtabITS), sampleinfoITS$Sample.ID),]
rownames(sampleinfoITS)=sampleinfoITS$Sample.ID
sampleinfoITS=sampleinfoITS[,-1]
OTU=otu_table(asvtabITS, taxa_are_rows = FALSE)
TAX=tax_table(as.matrix(TAXITS))
physeqmod = phyloseq(OTU, TAX, sample_data(sampleinfoITS))
#####save sequences and change names:
sequences=taxa_names(physeqmod)
ASVnames=sprintf("ASVfun%s",seq(1:ncol(physeqmod@otu_table)))
sequencesnames=rbind(ASVnames, sequences)
#saveRDS(sequencesnames, "/Users/gabri/OneDrive - University of Arizona/dust_project/data/soil/sequences_fungi.RDS")
taxa_names(physeqmod)=ASVnames






##### Remove playa point (point number 26).
##### Willcox playa is a dry lake in southern Arizona. It was sampled but given the high difference with the other sampled environments,
##### total lack of vegetation. Extremely high salinity, lack of organic carbon, very High pH
##### it will not be used since it is an outlier. 
physeqmod=prune_samples(sample_data(physeqmod)$point!=26, physeqmod)
physeqmod <- prune_taxa(taxa_sums(physeqmod) > 0, physeqmod) # Eliminate rows with only zeros











soil<- prune_samples(sample_data(physeqmod)$Type == "Soil", physeqmod)
###insert real metadata table
metadata_soil=read.csv("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/metadata_soil.csv", row.names = 1)
metadata_soil=sample_data(metadata_soil)
soil <- merge_phyloseq(soil, metadata_soil)
soil@sam_data=soil@sam_data[,-c(2,3,4)]
soil = prune_samples(sample_names(soil)!="DP29-02-S", soil) ###bad point
soil = prune_samples(sample_names(soil)!="DP03-01-S", soil) ###bad point
soil <- prune_taxa(taxa_sums(soil) > 0, soil)

soil_otu_fungi = data.frame(t(soil@otu_table))
soil_md_fungi =data.frame(soil@sam_data)
soil_taxa_fungi = data.frame(soil@tax_table)
#write.csv2(soil_otu_fungi, "/Users/gabri/Documents/otu_fungi.csv")
#write.csv2(soil_md_fungi, "/Users/gabri/Documents/md_fungi.csv")
#write.csv2(soil_taxa_fungi, "/Users/gabri/Documents/taxa_fungi.csv")




#in the dust part i need to add also those soil samples sampled in the dust places
dustpoints=c(1,10,11,18)
dust<- prune_samples(sample_data(physeqmod)$point %in% dustpoints, physeqmod)
metadata_dust=read.csv("/Users/gabri/Documents/GitHub/Dust_project/data/metadata/dust_met.csv", row.names = 1)
metadata_dust=sample_data(metadata_dust)
dust <- merge_phyloseq(dust, metadata_dust)
dust@sam_data=dust@sam_data[,-c(1,2,3)]


###summary of table
ASVtab_soil = data.frame(t(soil@otu_table))
nrow(ASVtab_soil)#4791
mean(colSums(ASVtab_soil!=0)) #number of sequences per sample 159.16
sd(colSums(ASVtab_soil!=0)) #SD 53.63
taxatab = data.frame(soil@tax_table)
ASVtab_soil$phylum = taxatab$Phylum[match(rownames(ASVtab_soil), rownames(taxatab))]
ASVtab_soil = aggregate(. ~ phylum ,data = ASVtab_soil,  FUN = sum)
rownames(ASVtab_soil) = ASVtab_soil[,1]
ASVtab_soil = ASVtab_soil[,-1]
rowsums_ASVtab_soil = rowSums(ASVtab_soil)
total_sum = sum(rowsums_ASVtab_soil)
fractions = rowsums_ASVtab_soil/total_sum
fractions = fractions[order(-fractions)]
fractions = fractions*100
fractions


###summary of table
ASVtab_dust = data.frame(t(dust@otu_table))
ASVtab_dust = ASVtab_dust[,-c(1,2,3,16,17,18,31,32,33,46,47,48)]
ASVtab_dust = ASVtab_dust[rowSums(ASVtab_dust[])>0,]
nrow(ASVtab_dust)#2515
mean(colSums(ASVtab_dust!=0)) #number of sequences per sample 181.19
sd(colSums(ASVtab_dust!=0)) #SD 70.8
taxatab = data.frame(dust@tax_table)
ASVtab_dust$phylum = taxatab$Phylum[match(rownames(ASVtab_dust), rownames(taxatab))]
ASVtab_dust = aggregate(. ~ phylum ,data = ASVtab_dust,  FUN = sum)
rownames(ASVtab_dust) = ASVtab_dust[,1]
ASVtab_dust = ASVtab_dust[,-1]
rowsums_ASVtab_dust = rowSums(ASVtab_dust)
total_sum = sum(rowsums_ASVtab_dust)
fractions = rowsums_ASVtab_dust/total_sum
fractions = fractions[order(-fractions)]
fractions = fractions*100
fractions





soil <- prune_taxa(taxa_sums(soil) > 0, soil)
saveRDS(soil, "/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soilITS_ASVnames.RDS")
dust <- prune_taxa(taxa_sums(dust) > 0, dust)
saveRDS(dust, "/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")

