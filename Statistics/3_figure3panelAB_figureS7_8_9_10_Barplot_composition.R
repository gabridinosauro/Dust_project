####barplot Composition
library(vegan)
library(phyloseq)
library(dplyr)# alternatively, this also loads %>%
library(ggplot2)
library(ggpubr)









####barplot16s soil####
soil=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soil16sASVnames.RDS")
taxatab = data.frame(soil@tax_table)
taxatab[is.na(taxatab)] <- 'unknown'
tax_table(soil) = as.matrix(taxatab)


dat.aglo = tax_glom(soil, taxrank = "Phylum")
dat.dataframe = psmelt(dat.aglo)
dat.dataframe$pH = as.numeric(dat.dataframe$pH )
dat.agr = aggregate(Abundance ~ point + Phylum, data=dat.dataframe, FUN=sum)
ph.aggr = aggregate(pH ~ point, data=dat.dataframe, FUN=mean)
dat.agr$Phylum= as.character(dat.agr$Phylum)
###keep 10 most abundant phyla
ten<-aggregate(Abundance~Phylum,data= dat.agr, FUN=sum)
ten<-ten[order(ten$Abundance, decreasing=TRUE),]
ten<-as.vector(ten$Phylum[1:10])
`%nin%` = Negate(`%in%`) # create a negation, for not included (there is not in R)
for ( i  in 1:length(dat.agr$Phylum)) #substitute those not in the top 10
  if(dat.agr$Phylum[i] %nin% ten) dat.agr$Phylum[i]<-"others"
dat.agr$pH =  ph.aggr[match(dat.agr$point, ph.aggr$point),2]
dat.agr$pH = round(dat.agr$pH, digits = 3)
names <- data.frame (points = c(1:25, 27:30),
                  renamed = c(1:29))
dat.agr$point = names[match(dat.agr$point, names$point),2]
dat.agr$point = as.factor(dat.agr$point)
dat.agr$point <- reorder(dat.agr$point, dat.agr$pH)

a<-ggplot(dat.agr, aes(x=point, y=Abundance, fill=Phylum)) + geom_bar(stat="identity", position = "fill") +
  theme_classic() +
  scale_fill_brewer(palette="Paired")  + theme(axis.text.x = element_text(angle =45, hjust = 1))  +
  xlab("Sampling location") + ylab("Reative abundance")
a
#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/3rd_draft/plots/SI/barplot_bacteria_soil.pdf", device = "pdf", height = 6, width = 10 )















####barplotITS soil####
soil=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_soilITS_ASVnames.RDS")
taxatab = data.frame(soil@tax_table)
taxatab[is.na(taxatab)] <- 'unknown'
tax_table(soil) = as.matrix(taxatab)

dat.aglo = tax_glom(soil, taxrank = "Phylum")
dat.dataframe = psmelt(dat.aglo)
dat.dataframe$pH = as.numeric(dat.dataframe$pH )
dat.agr = aggregate(Abundance ~ point + Phylum, data=dat.dataframe, FUN=sum)
ph.aggr = aggregate(pH ~ point, data=dat.dataframe, FUN=mean)
dat.agr$Phylum= as.character(dat.agr$Phylum)
###keep 10 most abundant phyla
ten<-aggregate(Abundance~Phylum,data= dat.agr, FUN=sum)
ten<-ten[order(ten$Abundance, decreasing=TRUE),]
ten<-as.vector(ten$Phylum[1:10])
`%nin%` = Negate(`%in%`) # create a negation, for not included (there is not in R)
for ( i  in 1:length(dat.agr$Phylum)) #substitute those not in the top 10
  if(dat.agr$Phylum[i] %nin% ten) dat.agr$Phylum[i]<-"others"
dat.agr$pH =  ph.aggr[match(dat.agr$point, ph.aggr$point),2]
dat.agr$pH = round(dat.agr$pH, digits = 3)
names <- data.frame (points = c(1:25, 27:30),
                     renamed = c(1:29))
dat.agr$point = names[match(dat.agr$point, names$point),2]
dat.agr$point = as.factor(dat.agr$point)
dat.agr$point <- reorder(dat.agr$point, dat.agr$pH)

a<-ggplot(dat.agr, aes(x=point, y=Abundance, fill=Phylum)) + geom_bar(stat="identity", position = "fill") +
  theme_classic() +
  scale_fill_brewer(palette="Paired")  + theme(axis.text.x = element_text(angle =45, hjust = 1))  +
  ylab("Relative abundance")
a











####barplotITS soi class####
dat.aglo = tax_glom(soil, taxrank = "Class")
dat.dataframe = psmelt(dat.aglo)
dat.dataframe$pH = as.numeric(dat.dataframe$pH )
dat.agr = aggregate(Abundance ~ point + Class, data=dat.dataframe, FUN=sum)
ph.aggr = aggregate(pH ~ point, data=dat.dataframe, FUN=mean)
dat.agr$Class= as.character(dat.agr$Class)
###keep 10 most abundant phyla
ten<-aggregate(Abundance~Class,data= dat.agr, FUN=sum)
ten<-ten[order(ten$Abundance, decreasing=TRUE),]
ten<-as.vector(ten$Class[1:10])
`%nin%` = Negate(`%in%`) # create a negation, for not included (there is not in R)
for ( i  in 1:length(dat.agr$Class)) #substitute those not in the top 10
  if(dat.agr$Class[i] %nin% ten) dat.agr$Class[i]<-"others"
dat.agr$pH =  ph.aggr[match(dat.agr$point, ph.aggr$point),2]
dat.agr$pH = round(dat.agr$pH, digits = 3)
names <- data.frame (points = c(1:25, 27:30),
                     renamed = c(1:29))
dat.agr$point = names[match(dat.agr$point, names$point),2]
dat.agr$point = as.factor(dat.agr$point)
dat.agr$point <- reorder(dat.agr$point, dat.agr$pH)

b<-ggplot(dat.agr, aes(x=point, y=Abundance, fill=Class)) + geom_bar(stat="identity", position = "fill") +
  theme_classic() +
  scale_fill_brewer(palette="Paired")  + theme(axis.text.x = element_text(angle =45, hjust = 1))  +
  ylab("Relative abundance")
b
ggarrange(a,b, ncol = 1)
#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/3rd_draft/plots/SI/barplot_fungi_phyla_soil.pdf", device = "pdf", height = 12, width = 10 )























####barplot16s dust####
dust=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
taxatab = data.frame(dust@tax_table)
taxatab[is.na(taxatab)] <- 'unknown'
tax_table(dust) = as.matrix(taxatab)
dat.aglo = tax_glom(dust, taxrank = "Phylum")
#dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.aglo)
dat.agr = aggregate(Abundance~Dust_Type+Phylum+Point_type, data=dat.dataframe, FUN=sum)
dat.agr$Phylum= as.character(dat.agr$Phylum)
###keep 10 most abundant phyla
ten<-aggregate(Abundance~Phylum,data= dat.agr, FUN=sum)
ten<-ten[order(ten$Abundance, decreasing=TRUE),]
ten<-as.vector(ten$Phylum[1:10])
`%nin%` = Negate(`%in%`) # create a negation, for not included (there is not in R)
for ( i  in 1:length(dat.agr$Phylum)) #substitute those not in the top 10
  if(dat.agr$Phylum[i] %nin% ten) dat.agr$Phylum[i]<-"others"
dat.agr$Dust_Type <- factor(dat.agr$Dust_Type ,levels = c("Soil","500","75", "25", "Dust_Gen"))
dat.agr$Point_type <- factor(dat.agr$Point_type)
levels(dat.agr$Point_type)  <- c("DS1", "DS2", "PP", "SCR")


a<-ggplot(dat.agr, aes(x=Dust_Type, y=Abundance, fill=Phylum)) + geom_bar(stat="identity", position = "fill") +
  theme_classic() +
  scale_fill_brewer(palette="Paired")  + theme(axis.text.x = element_text(angle =45, hjust = 1))  + facet_wrap("Point_type", nrow = 1,  scales="free")+
  scale_x_discrete(labels=c("Soil", "Coarse LP", "Medium LP", "Fine LP", "WT")) + ylab("Relative Abundance") 
a



##########Condense in one point
dust=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dust16sASVnames.RDS")
taxatab = data.frame(dust@tax_table)
taxatab[is.na(taxatab)] <- 'unknown'
tax_table(dust) = as.matrix(taxatab)

dat.aglo = tax_glom(dust, taxrank = "Phylum")
#dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.aglo)
dat.agr = aggregate(Abundance~Dust_Type+Phylum, data=dat.dataframe, FUN=sum)
dat.agr$Phylum= as.character(dat.agr$Phylum)
###keep 10 most abundant phyla
ten<-aggregate(Abundance~Phylum,data= dat.agr, FUN=sum)
ten<-ten[order(ten$Abundance, decreasing=TRUE),]
ten<-as.vector(ten$Phylum[1:10])
`%nin%` = Negate(`%in%`) # create a negation, for not included (there is not in R)
for ( i  in 1:length(dat.agr$Phylum)) #substitute those not in the top 10
  if(dat.agr$Phylum[i] %nin% ten) dat.agr$Phylum[i]<-"others"
dat.agr$Dust_Type <- factor(dat.agr$Dust_Type ,levels = c("Soil","500","75", "25", "Dust_Gen"))


a_bac<-ggplot(dat.agr, aes(x=Dust_Type, y=Abundance, fill=Phylum)) + geom_bar(stat="identity", position = "fill") +
  theme_classic() +
  scale_fill_brewer(palette="Paired")  + theme(axis.text.x = element_text(angle =45, hjust = 1)) +
  scale_x_discrete(labels=c("Soil", "Coarse LP", "Medium LP", "Fine LP", "WT")) + ylab("Relative Abundance") + xlab(NULL) +
  ggtitle("Bacteria/Archaea")














####barplotITS dust####
dust=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
taxatab = data.frame(dust@tax_table)
taxatab[is.na(taxatab)] <- 'unknown'
tax_table(dust) = as.matrix(taxatab)

dat.aglo = tax_glom(dust, taxrank = "Phylum")
#dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.aglo)
dat.agr = aggregate(Abundance~Dust_Type+Phylum+Point_type, data=dat.dataframe, FUN=sum)
dat.agr$Phylum= as.character(dat.agr$Phylum)
###keep 10 most abundant phyla
ten<-aggregate(Abundance~Phylum,data= dat.agr, FUN=sum)
ten<-ten[order(ten$Abundance, decreasing=TRUE),]
ten<-as.vector(ten$Phylum[1:10])
`%nin%` = Negate(`%in%`) # create a negation, for not included (there is not in R)
for ( i  in 1:length(dat.agr$Phylum)) #substitute those not in the top 10
  if(dat.agr$Phylum[i] %nin% ten) dat.agr$Phylum[i]<-"others"
dat.agr$Dust_Type <- factor(dat.agr$Dust_Type ,levels = c("Soil","500","75", "25", "Dust_Gen"))
dat.agr$Point_type <- factor(dat.agr$Point_type)
levels(dat.agr$Point_type)  <- c("DS1", "DS2", "PP", "SCR")


a<-ggplot(dat.agr, aes(x=Dust_Type, y=Abundance, fill=Phylum)) + geom_bar(stat="identity", position = "fill") +
  theme_classic() +
  scale_fill_brewer(palette="Paired")  + theme(axis.text.x = element_text(angle =45, hjust = 1))  + facet_wrap("Point_type", nrow = 1,  scales="free")+
  scale_x_discrete(labels=c("Soil", "Coarse LP", "Medium LP", "Fine LP", "WT")) + ylab("Relative Abundance") + xlab("Sample type")
a

dust=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
taxatab = data.frame(dust@tax_table)
taxatab[is.na(taxatab)] <- 'unknown'
tax_table(dust) = as.matrix(taxatab)
dat.aglo = tax_glom(dust, taxrank = "Class")
#dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.aglo)
dat.agr = aggregate(Abundance~Dust_Type+Class+Point_type, data=dat.dataframe, FUN=sum)
dat.agr$Class= as.character(dat.agr$Class)
###keep 10 most abundant phyla
ten<-aggregate(Abundance~Class,data= dat.agr, FUN=sum)
ten<-ten[order(ten$Abundance, decreasing=TRUE),]
ten<-as.vector(ten$Class[1:10])
`%nin%` = Negate(`%in%`) # create a negation, for not included (there is not in R)
for ( i  in 1:length(dat.agr$Class)) #substitute those not in the top 10
  if(dat.agr$Class[i] %nin% ten) dat.agr$Class[i]<-"others"
dat.agr$Dust_Type <- factor(dat.agr$Dust_Type ,levels = c("Soil","500","75", "25", "Dust_Gen"))
dat.agr$Point_type <- factor(dat.agr$Point_type)
levels(dat.agr$Point_type)  <- c("DS1", "DS2", "PP", "SCR")


b<-ggplot(dat.agr, aes(x=Dust_Type, y=Abundance, fill=Class)) + geom_bar(stat="identity", position = "fill") +
  theme_classic() +
  scale_fill_brewer(palette="Paired")  + theme(axis.text.x = element_text(angle =45, hjust = 1))  + facet_wrap("Point_type", nrow = 1,  scales="free")+
  scale_x_discrete(labels=c("Soil", "Coarse LP", "Medium LP", "Fine LP", "WT")) + ylab("Relative Abundance") + xlab("Sample type")
b


ggarrange(a,b, ncol = 1, labels = "AUTO")
#ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/3rd_draft/plots/SI/barplot_fungi_phyla_dust.pdf", device = "pdf", height = 12, width = 10 )




########make them just one
dust=readRDS("/Users/gabri/Documents/GitHub/Dust_project/data/phyloseq_obects/phyloseq_dustITS_ASVnames.RDS")
taxatab = data.frame(dust@tax_table)
taxatab[is.na(taxatab)] <- 'unknown'
tax_table(dust) = as.matrix(taxatab)
dat.aglo = tax_glom(dust, taxrank = "Class")
#dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.aglo)
dat.agr = aggregate(Abundance~Dust_Type+Class, data=dat.dataframe, FUN=sum)
dat.agr$Class= as.character(dat.agr$Class)
###keep 10 most abundant phyla
ten<-aggregate(Abundance~Class,data= dat.agr, FUN=sum)
ten<-ten[order(ten$Abundance, decreasing=TRUE),]
ten<-as.vector(ten$Class[1:10])
`%nin%` = Negate(`%in%`) # create a negation, for not included (there is not in R)
for ( i  in 1:length(dat.agr$Class)) #substitute those not in the top 10
  if(dat.agr$Class[i] %nin% ten) dat.agr$Class[i]<-"others"
dat.agr$Dust_Type <- factor(dat.agr$Dust_Type ,levels = c("Soil","500","75", "25", "Dust_Gen"))



b<-ggplot(dat.agr, aes(x=Dust_Type, y=Abundance, fill=Class)) + geom_bar(stat="identity", position = "fill") +
  theme_classic() +
  scale_fill_brewer(palette="Paired")  + theme(axis.text.x = element_text(angle =45, hjust = 1)) +
  scale_x_discrete(labels=c("Soil", "Coarse SS", "Medium SS", "Fine SS", "WT")) + ylab("Relative Abundance")  + xlab(NULL)+
  ggtitle("Fungi")
b

ggarrange(a_bac,b, ncol = 1, labels = "AUTO")
ggsave("/Users/gabri/OneDrive - University of Arizona/dust_project/3rd_draft/plots/figure_3/barplots.pdf", device = "pdf", height = 7, width = 3.5 )




