library(dada2)
getwd()
setwd("/groups/barberan/Gabri/dp/16s/Fastqfiles16S/demuinvpaltes")
getwd()
path <- "/groups/barberan/Gabri/dp/16s/Fastqfiles16S/demuinvpaltes"
list.files(path, pattern = "*fastq*")
fnFs <- sort(list.files(path, pattern="_R1", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2", full.names = TRUE))
sample.names0 <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
sample.names = gsub(".fastq.gz", "", sample.names0)
plotQualityProfile(fnFs[11:20]) # forward
plotQualityProfile(fnRs[12:20]) # reverse  
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# Trim forwards at 140bp and reverses at 140bp based on quality plots; also trim first 10 bp as recommended
# After this completes, check "out" to see if filtering was too stringent
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft = 10) 
# Monitor number of sequences
out = as.data.frame(out)
out$Retained = out$reads.out/out$reads.in*100
head(out)
out$filtFs<-filtFs
out$filtRs<-filtRs
out1<-out[which(out$reads.out!=0), ] 
filtFs<-out1$filtFs
filtRs<-out1$filtRs
out<-out[,-c(4,5)]
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#save.image("aftererror.RData")

# Reduce computation time by dereplicating sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
sample.names <- sapply(strsplit(basename(filtFs), "_F"), `[`, 1)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
save.image("pre_error.RData")
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "Retained %", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
taxa <- assignTaxonomy(seqtab.nochim, "/home/u1/schiro/silva_nr_v132_train_set.fa.gz?download=1", multithread=TRUE)
taxa <- addSpecies(taxa, "/home/u1/schiro/silva_species_assignment_v132.fa.gz?download=1")
saveRDS(taxa, file = "TAX16sinvplates.rds")#taxatable
saveRDS(seqtab.nochim, file = "OTU16sinvplates.RDS")#OTu table
write.table(track, "sequence_pipeline_stats.txt", sep = "\t", quote = F)
#save.image("finished.RData")
#load("finished.RData")
