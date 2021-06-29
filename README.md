# Dust project

Hello,
This is a step by step process to repeat the data analysis done in the paper:
"Linking dust dispersal and spatial distribution of microbes in a complex semi-arid region"

I will first have to download the data.
There are two options:
1)You can download the whole raw data from this link:
2)Or you can dowloand the data combined phyloseq objects (see links below).


There are three steps: 
1)The first one is the bioinformatic procedure that, from demultiplexed data brings you to two ASV count tables one Bacteria/Archaea and one Fungi (Dada2 pipeline). Taxa tables will be also produced.

2)The second step brings you from an ASV count table, for both Bacteria/Archaea communities to a  phyloseq objects where ASV count tables, taxonomy and metadata of each condition is going to be used.

For convenience, the data will be divided in four data sets.

(1) Spatial dataset for Bacteria/Archaea, this includes all the soil samples analyzed,
(2) Spatial dataset for Fungi,
(3) Dust dataset for Bacteria/Archaea,  this includes soil samples and dust samples from the locations PP, SCR, DS1 and DS2,
(4) Dust dataset for Bacteria/Archaea,  this includes soil samples and dust samples from the locations PP, SCR, DS1 and DS2.

3)The third step includes the statistical analyses done.

