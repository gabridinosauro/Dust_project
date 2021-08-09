# Dust project

Hello,
This is a step by step process to repeat the data analysis done in the paper:
"Linking dust dispersal and spatial distribution of microbes in a complex semi-arid region"

first you have to download the data.
There are two options:

1. You can download the whole raw data from this link INSERT LINK: (highly time and memory consuming (GB of data).
2. Or you can download the combined phyloseq objects (see links below) (quick and fast).


There are three steps: 
1. Only from raw data (point one above) the bioinformatic procedure that, from demultiplexed data brings you to two ASV count tables one Bacteria/Archaea and one Fungi (Dada2 pipeline). Taxa tables will be also produced. Scripts are inside the folder "bioinformatics". End data will be saved as RDS file, to be easily imported into R. The final results are all inside the {"end_dada2" folder). These files are necessary only for the following step.

2. The second step brings you from the outputs of dada2 for both Bacteria/Archaea communities to phyloseq objects where ASV count tables, taxonomy and metadata of each condition are cleant and organized in the two different datasets (dust and soil).
The output will therefore create the four data sets used for statistical analysis: <br />

      a)Spatial dataset for Bacteria/Archaea, this includes all the soil samples analyzed,<br />
      b)Spatial dataset for Fungi,<br />
      c)Dust dataset for Bacteria/Archaea,  this includes soil samples and dust samples from the locations PP, SCR, DS1 and DS2,<br />
      d)Dust dataset for Fungi,  this includes soil samples and dust samples from the locations PP, SCR, DS1 and DS2.<br />

3. The third step includes the statistical analyses done, consisting in various scripts from 1 to 15, in the stats folder, recreating step by step the analysis done in the paper. (script number 2 "2_geographic_calculations.R" will produce the file "geo_dist_data.RDS" and "geo_dist_dataITS.RDS" necessary for some of the scripts;
scripts "11_DEseq2_bacteria_and_heatmaps.R " and "12_DEseq2_ITS_and_heatmaps.R" will create all the files within the folder "data/Deseq2", also necessary for some of the scripts.

