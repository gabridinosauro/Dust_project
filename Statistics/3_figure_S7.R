library(ggplot2)
library(vegan)
library(phyloseq)
library(pheatmap)
library(viridis)
library(DESeq2)
library(pheatmap)
library(metagMisc)
library(rcompanion)
library(car)
library(multcompView)
library(lsmeans)
library(multcomp)
library(agricolae)
library(ggpubr)

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
sam_data$Dust_Type[sam_data$Dust_Type == "500"] <- "Coarse SS"
sam_data$Dust_Type[sam_data$Dust_Type == "75"] <- "Medium SS"
sam_data$Dust_Type[sam_data$Dust_Type == "25"] <- "Fine SS"
sam_data$Dust_Type<- factor(sam_data$Dust_Type,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))







###calculate community indexes across sample types
summary(rowSums(t(asvtab1)))
sam_data$Richness.rar<-specnumber(rrarefy(t(asvtab1), 14500))
sam_data$Shannon.rar = vegan::diversity(rrarefy(t(asvtab1), 14500))
sam_data$Evenness.rar = sam_data$Shannon.rar/log(sam_data$Richness.rar)





####statistical test richness
model = lm(Richness.rar ~ Dust_Type + Point_type, data = sam_data) 
anova(lm(Richness.rar ~ Dust_Type + Point_type, data = sam_data))
model.av = aov(model)
tukey.test = HSD.test(model.av, trt = "Dust_Type")
tukey.test



####plot richness
richness_bac = ggplot(sam_data, aes(x=Dust_Type, y=Richness.rar, fill  = Dust_Type)) + 
  geom_point(position = "jitter", aes(color = Dust_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5)  + 
  theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","b","b","b","b" ), fun = max, vjust = 0) +ylab("Number of phylotypes") +xlab(NULL)+
  ggtitle("Richness")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
richness_bac





###statistical test diversity
model = lm(Shannon.rar ~ Dust_Type + Point_type, data = sam_data) 
anova(lm(Shannon.rar ~ Dust_Type + Point_type, data = sam_data))
model.av = aov(model)
tukey.test = HSD.test(model.av, trt = "Dust_Type")
tukey.test



###plot diversity
div_bac = ggplot(sam_data, aes(x=Dust_Type, y=Shannon.rar, fill  = Dust_Type)) + 
  geom_point(position = "jitter", aes(color = Dust_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5)  + 
  theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","a","a","a","a" ), fun = max, vjust = 0) +ylab("Shannon Index (H)") +xlab(NULL)+
  ggtitle("Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
div_bac







####statistical test evenness
#shapiro.test(sam_data$Evenness.rar)
model = lm(Evenness.rar ~ Dust_Type + Point_type, data = sam_data) 
anova(lm(Evenness.rar ~ Dust_Type + Point_type, data = sam_data))
model.av = aov(model)
par(mfrow=c(2,2))
#plot(model)
tukey.test = HSD.test(model.av, trt = "Dust_Type")
tukey.test

evenness_bac = ggplot(sam_data, aes(x=Dust_Type, y=Evenness.rar, fill  = Dust_Type)) + 
  geom_point(position = "jitter", aes(color = Dust_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5)  + 
  theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("a","b","b","b","b" ), fun = max, vjust = 0) +ylab("Pielou's evenness  (J)") +xlab(NULL)+
  ggtitle("Evenness") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
evenness_bac


a  = ggarrange(richness_bac, div_bac, evenness_bac, nrow = 1, labels = "AUTO")
a = annotate_figure(a, top = text_grob("Bacteria/Archaea"))
a




####Fungi######

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
sam_data$Dust_Type[sam_data$Dust_Type == "500"] <- "Coarse SS"
sam_data$Dust_Type[sam_data$Dust_Type == "75"] <- "Medium SS"
sam_data$Dust_Type[sam_data$Dust_Type == "25"] <- "Fine SS"
sam_data$Dust_Type<- factor(sam_data$Dust_Type,levels = c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT"))


summary(rowSums(t(asvtab1)))
sam_data$Richness.rar<-specnumber(rrarefy(t(asvtab1), 15500))
sam_data$Shannon.rar = vegan::diversity(rrarefy(t(asvtab1), 15500))
sam_data$Evenness.rar = sam_data$Shannon.rar/log(sam_data$Richness.rar)






model = lm(Richness.rar ~ Dust_Type + Point_type, data = sam_data) 
anova(lm(Richness.rar ~ Dust_Type + Point_type, data = sam_data))
model.av = aov(model)
tukey.test = HSD.test(model.av, trt = "Dust_Type")
tukey.test

richness_fung = ggplot(sam_data, aes(x=Dust_Type, y=Richness.rar, fill  = Dust_Type)) + 
  geom_point(position = "jitter", aes(color = Dust_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5)  + 
  theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("b","a","a","a","ab" ), fun = max, vjust = 0) +ylab("Number of phylotypes") +xlab(NULL)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Richness")
richness_fung


model = lm(Shannon.rar ~ Dust_Type + Point_type, data = sam_data) 
anova(lm(Shannon.rar ~ Dust_Type + Point_type, data = sam_data))
model.av = aov(model)
tukey.test = HSD.test(model.av, trt = "Dust_Type")
tukey.test

div_fun = ggplot(sam_data, aes(x=Dust_Type, y=Shannon.rar, fill  = Dust_Type)) + 
  geom_point(position = "jitter", aes(color = Dust_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5)  + 
  theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("ab","a","ab","ab","b" ), fun = max, vjust = 0) +ylab("Shannon Index (H)") +xlab(NULL)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Diversity")
div_fun





model = lm(Evenness.rar ~ Dust_Type + Point_type, data = sam_data) 
anova(lm(Evenness.rar ~ Dust_Type + Point_type, data = sam_data))
model.av = aov(model)
tukey.test = HSD.test(model.av, trt = "Dust_Type")
tukey.test

evenness_fun = ggplot(sam_data, aes(x=Dust_Type, y=Evenness.rar, fill  = Dust_Type)) + 
  geom_point(position = "jitter", aes(color = Dust_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5)  + 
  theme_classic() + theme(legend.position = "none") +
  stat_summary(geom = 'text', label = c("ab","a","ab","ab","b" ), fun = max, vjust = 0) +ylab("Pielou's evenness  (J)") +xlab(NULL)+
  ggtitle("Evenness") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
evenness_fun


b  = ggarrange(richness_fung, div_fun, evenness_fun, nrow = 1, labels = c("D","E","F"))
b = annotate_figure(b, top = text_grob("Fungi"))

ggarrange(a,b, nrow = 2)




