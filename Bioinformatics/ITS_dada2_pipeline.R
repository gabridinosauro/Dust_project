library(dada2)
library(ShortRead)
library(Biostrings)
path <- "/groups/barberan/Gabri/dp/ITS/demultiplexedITS"
list.files(path, pattern = "*fastq*")
# Assign files to objects
fnFs <- sort(list.files(path, pattern="_R1", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2", full.names = TRUE))
#Identify primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  # forward primer sequence
REV <- "GCTGCGTTCTTCATCGATGC"  # reverse primer sequence

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

FWD.orients
REV.orients
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[100]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[100]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[100]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[100]]))


sample.names0 <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
sample.names = gsub(".fastq.gz", "", sample.names0)

# Inspect quality of sequences in files
plotQualityProfile(fnFs[35:40]) # forward
plotQualityProfile(fnRs[35:40]) # reverse  
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), trimRight = c(0, 20), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)
out
#Learn the error rates
filtFs<- filtFs[-c(9,10)]
errF <- learnErrors(filtFs[5:18], multithread = TRUE)
filtRs<- filtRs[-c(9,10)]
errR <- learnErrors(filtRs[5:18], multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names[-c(9,10)]
names(derepRs) <- sample.names[-c(9,10)]
#save.image("/groups/barberan/Gabri/dp/ITS/workspaces/ITS.RData")
#load("/groups/barberan/Gabri/dp/ITS/workspaces/ITS.RData")
#Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# See distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim)))

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
head(track)
# Download
system("wget -N https://files.plutof.ut.ee/public/orig/E7/28/E728E2CAB797C90A01CD271118F574B8B7D0DAEAB7E81193EB89A2AC769A0896.gz")
system("gunzip E728E2CAB797C90A01CD271118F574B8B7D0DAEAB7E81193EB89A2AC769A0896.gz")
#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "sh_general_release_dynamic_04.02.2020.fasta", multithread=TRUE, tryRC=T)
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
View(taxa.print)
write.table(taxa, "taxa_its.txt", sep="\t", quote=F, row.names=T)
write.table(seqtab.nochim, "asv_its.txt", sep="\t", quote=F, row.names=T)
write.table(track, "sequence_stats_its.txt", sep="\t", quote=F, row.names=T)
#save.image("/groups/barberan/Gabri/dp/ITS/workspaces/ITS.RData")
