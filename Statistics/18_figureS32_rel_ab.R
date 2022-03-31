library(phyloseq)


#### Plot the curve I want to plot###
#########LOAD AND PREPARE THE DATA##################
asvtab<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
ASVtab_dust <- prune_taxa(taxa_sums(asvtab) > 0, asvtab)#elimiante_zero_counts
rank_names(ASVtab_dust)
sam_data = data.frame(ASVtab_dust@sam_data)
sam_data = sam_data[,-1] #eliminate first column (useless)
sam_data$sample <- paste(sam_data$Point_type,sam_data$Repetition, sep = "_")






#####insert new asvtab into the phyloseq object
ASVtab_dust@sam_data =sample_data(sam_data) 
ASVtab_dust@sam_data$Repetition=factor(ASVtab_dust@sam_data$Repetition)
rm(asvtab)
asvtab1 = data.frame(t(ASVtab_dust@otu_table))
sam_data$Dust_Type[sam_data$Dust_Type == "Dust_Gen"] <- "WT"
sam_data$Dust_Type[sam_data$Dust_Type == "500"] <- "Coarse LP"
sam_data$Dust_Type[sam_data$Dust_Type == "75"] <- "Medium LP"
sam_data$Dust_Type[sam_data$Dust_Type == "25"] <- "Fine LP"
sam_data$Dust_Type<- factor(sam_data$Dust_Type,levels = c("Soil", "Coarse LP", "Medium LP", "Fine LP", "WT"))










asv_tab_rar = t(rrarefy(t(asvtab1), 14500))
asvtab1_rel = decostand(asv_tab_rar, method = "total", MARGIN = 2)


##frequency, test different thresholds, quick and dirty way)
df = data.frame()
for (i in 5:9)
{
  a = exp(-i)
  rare = apply(asvtab1_rel,2,function(x){sum(x < a & x>0)})
  total = apply(asvtab1_rel,2,function(x){sum( x>0)})
  freq = rare/total
  partial = data.frame(freq = freq, names = names(freq), numer = a)
  df = rbind(df, partial)
  
}


df$names = gsub("\\.", "-", df$names)
#$numer = log(df$numer)
df$numer = factor(df$numer, levels =unique(df$numer) )
df$Sample_type = sam_data$Dust_Type[match(df$names ,rownames(sam_data) )]
p<-ggplot(df, aes(x=numer, y=freq, fill=Sample_type)) +
  geom_boxplot() +theme_classic() +xlab("Cutoff relative abundance") +ylab("Relative Frequency") +
  scale_x_discrete(labels = c("e^-4","e^-5","e^-6","e^-7","e^-8","e^-7","e^-8","e^-9","e^-10"))+
  ggtitle("Archaea/Bacteria") 

p






###########FUNGI##############
#########LOAD AND PREPARE THE DATA##################
ASVtab_dust_fungi<-readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
sam_data = data.frame(ASVtab_dust_fungi@sam_data)
sam_data = sam_data[,-1] #eliminate first column (useless)
sam_data$sample <- paste(sam_data$Point_type,sam_data$Repetition, sep = "_")


#####insert new asvtab into the phyloseq object
ASVtab_dust_fungi@sam_data =sample_data(sam_data) 
ASVtab_dust_fungi@sam_data$Repetition=factor(ASVtab_dust_fungi@sam_data$Repetition)
rm(asvtab)
asvtab1 = data.frame(t(ASVtab_dust_fungi@otu_table))
sam_data$Dust_Type[sam_data$Dust_Type == "Dust_Gen"] <- "WT"
sam_data$Dust_Type[sam_data$Dust_Type == "500"] <- "Coarse LP"
sam_data$Dust_Type[sam_data$Dust_Type == "75"] <- "Medium LP"
sam_data$Dust_Type[sam_data$Dust_Type == "25"] <- "Fine LP"
sam_data$Dust_Type<- factor(sam_data$Dust_Type,levels = c("Soil", "Coarse LP", "Medium LP", "Fine LP", "WT"))



asv_tab_rar = t(rrarefy(t(asvtab1), 15500))
asvtab1_rel = decostand(asv_tab_rar, method = "total", MARGIN = 2)


##frequency, test different thresholds, quick and dirty way)
df = data.frame()
for (i in 5:9)
{
  a = exp(-i)
  rare = apply(asvtab1_rel,2,function(x){sum(x < a & x>0)})
  total = apply(asvtab1_rel,2,function(x){sum( x>0)})
  freq = rare/total
  partial = data.frame(freq = freq, names = names(freq), numer = a)
  df = rbind(df, partial)
  
}


df$names = gsub("\\.", "-", df$names)
df$numer = factor(df$numer, levels =unique(df$numer) )
df$Sample_type = sam_data$Dust_Type[match(df$names ,rownames(sam_data) )]
fun<-ggplot(df, aes(x=numer, y=freq, fill=Sample_type)) +
  geom_boxplot() +theme_classic() +xlab("Cutoff relative abundance") +ylab("Relative Frequency") +
  scale_x_discrete(labels = c("e^-4","e^-5","e^-6","e^-7","e^-8","e^-6","e^-7","e^-8","e^-9","e^-10"))+
  ggtitle("Fungi")

fun

ggarrange(p,fun, common.legend = TRUE)


