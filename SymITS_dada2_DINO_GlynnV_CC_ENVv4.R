library("knitr")
library(ggplot2)
library(gridExtra)
library(dada2)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(tidyverse)

knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
set.seed(100)


## -----------------------------------------------------------------------------
#Path will change for each run
path <- "/home/vmglynn/projects/def-barrett/vmglynn/internship2020/FULL_DINO/Seq"
list.files(path)

# Removed 19 files which only have F read on NCBI: SRR3193735 SRR3193738 SRR3193742 SRR3193743 SRR7079975 SRR7079976 SRR7079977 SRR7080011 SRR7080013 SRR7080014 SRR7080062 SRR7080063 SRR7080064 SRR7080078 SRR7080084 SRR7080085 SRR7080092 SRR7080093 SRR7080095



## -----------------------------------------------------------------------------

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)



## -----------------------------------------------------------------------------
plotQualityProfile(fnFs[1:2])
#Trim at 200, as ITS-2 not shorter than ~200 bp or longer than ~300 bp length (Hume et al. 2018)  and reads are exceptionally good quality. DADA2 is quite stringent, so better keep more than less reads at this step.

#Yet, as the ITS2 primer sets are not the same length, will need to trim F and R reads differently to a fixed length.  


## -----------------------------------------------------------------------------
plotQualityProfile(fnRs[1:2])
#Trim at 200, as ITS-2 not shorter than ~200 bp or longer than ~300 bp length (Hume et al. 2018) and reads are exceptionally good quality. DADA2 is quite stringent, so better keep more than less reads at this step.


## -----------------------------------------------------------------------------
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }


## ----filter and trim----------------------------------------------------------

# Place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered1", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered2", paste0(sample.names, "_R_filt.fastq.gz"))

any(duplicated(c(fnFs, fnRs)))
any(duplicated(c(filtFs, filtRs)))

length(fnFs)
length(fnRs)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=22, trimRight=28, truncLen=c(200,200),
             maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
             compress=TRUE, multithread=TRUE)
head(out)


## ----dereplication step-------------------------------------------------------

## Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence. Dereplication substantially reduces computation time by eliminating redundant comparisons.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

## verbose=TRUE, results in the service generating more output (will show you both WARNING and INFO log levels), normally you will only see WARNING or higher (ERROR for example).

# Name the derep-class objects by the sample names

names(derepFs) <- sample.names
names(derepRs) <- sample.names




## ----estimate error-----------------------------------------------------------

#Estimate error from first 40 samples
## The learnErrors method learns from a parametric error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution. 
## As in many machine-learning problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).

errF <- learnErrors(derepFs, multithread=TRUE)
errR <- learnErrors(derepRs, multithread=TRUE)



## ----plot errors--------------------------------------------------------------

# Plot estimated error rates, sanity check

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

## nominalQ, if TRUE, plot the expected error rates (red line) if quality scores exactly matched their nominal definition: Q = -10 log10(p_err).

## The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. 



## ----pseudopooling approach---------------------------------------------------
# Core sample inference algorithm

## Execution time

system.time(dd.pseudo <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo", verbose=0))

system.time(dd.pseudo <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo", verbose=0))

### Notice that in terms of processing, time from least to greatest: no pool, pseudo-pool, true pool. 

## Compare execution time with true pool and no pool 

### No pool
system.time(dd.pseudo <- dada(derepFs, err=errF, multithread=TRUE, pool=FALSE, verbose=0))

### True pool
system.time(dd.pseudo <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE, verbose=0))

# Sample inference with pseudo-pooling

dadaFs <- dada(derepFs, errF, pool=TRUE, multithread = TRUE)
dadaRs <- dada(derepRs, errR, pool=TRUE, multithread = TRUE)

dadaFs[[1]]
dadaRs[[1]]

## Pseudo-pooling is where samples are processed independently after sharing information between samples, approximating pooled sample inference in linear time.

## To pool or not to pool? Rarity biological-relevance vs PCR-artifact.

## OMEGA_A:  The key sensitivity parameters, controls the p-value threshold at which to call new ASVs.

## OMEGA_C: The error-correction threshold. One alternative is to turn off error-correction.
 


## ----merge paired reads-------------------------------------------------------
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Inspect merger data.frame. from first sample

head(mergers[[1]])
mergers[[1]]

## Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region.

## Most of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads?



## ----sequence table-----------------------------------------------------------

#Construct Sequence Table 
seqtab <- makeSequenceTable(mergers)

#Explore size of seqs
dim(seqtab)

## Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

hist(nchar(getSequences(seqtab)))

hist(nchar(getSequences(seqtab)), 
     main="Length distribution of seqs", 
     xlab="ITS-2 Size", 
     xlim=c(200,450),
     las=1, 
     breaks=50)


## ----remove chimeras----------------------------------------------------------

## Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

## Sequences that are much longer or shorter than expected may be the result of non-specific priming. You can remove non-target-length sequences with base R manipulations of the sequence table (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)]). This is analogous to “cutting a band” in-silico to get amplicons of the targeted length.



## ----seq length filter--------------------------------------------------------

# Sequence table length filter

MINLEN <- 150
MAXLEN <- 350

# MIN/MAXLEN can change based on histogram plot 

seqlens <- nchar(getSequences(seqtab.nochim)) 

st_len_filt <- seqtab.nochim[,seqlens >= MINLEN & seqlens <= MAXLEN]

## st_len_filt <- seqtab.nochim(seqlens = 150) 
                             
## st_len_filt <- seqtab.nochim(seqlens = 350) 

hist(nchar(getSequences(st_len_filt)), 
     main="Length Distribution of seqs", 
     xlab="ITS-2 Size", 
     xlim=c(150,350),
     las=1, 
     breaks=50)


## ----seq abundance filter-----------------------------------------------------

# Sequence table abundance filter

MINABUND <- 1 # set to 1 perhaps, can change
abundances <- colSums(st_len_filt)
st_len_abu_filt <- st_len_filt[,abundances >= MINABUND]

hist(nchar(getSequences(st_len_abu_filt)), 
     main="Abundance distribution of seqs", 
     xlab="Abundance", 
     xlim=c(150, 400),
     las=1, 
     breaks=50)


## ----further exploration------------------------------------------------------

seqtab.nochim2 <- as_tibble(seqtab.nochim)
glimpse(seqtab.nochim2)
str(seqtab.nochim2)

#Return sequences with col 1 with no heading but number and column 2 heading = x

seqs.nochim <- getSequences(seqtab.nochim)
head(seqs.nochim)
str(seqs.nochim)

## Compactly display the internal structure of an R object, a diagnostic function and an alternative to summary (and to some extent, dput). Ideally, only one line for each ‘basic’ structure is displayed. It is especially well suited to compactly display the (abbreviated) contents of (possibly nested) lists; str(object, …) = any R object about which you want to have some information.

## head() = return first part of vector, matrix, table, df, or function. 

## remove: colnames(seqs.nochim) <- c("sequence")
### Received "Error in `colnames<-`(`*tmp*`, value = "sequence") : attempt to set 'colnames' on an object with less than two dimensions"

# Write seq table

write.csv(seqtab.nochim, file = "seqtab.nochim.csv")
write.csv(st_len_filt, file = "st_len_filt.csv")
write.csv(st_len_abu_filt, file = "st_len_abu_filt.csv")


## ----compare filter steps 1---------------------------------------------------
hist(rowSums(seqtab.nochim))
hist(rowSums(st_len_filt))
hist(rowSums(st_len_abu_filt))

row_seqtab = tibble(rowSums(seqtab.nochim))
row_slaf = tibble(rowSums(st_len_abu_filt))

qplot(row_seqtab$`rowSums(seqtab.nochim)`, geom="histogram", binwidth=5000)
qplot(row_slaf$`rowSums(st_len_abu_filt)`, geom="histogram", binwidth=5000)

hist(colSums(st_len_filt))
hist(colSums(st_len_abu_filt))
## REMINDER: abundances <- colSums(st_len_filt) | st_len_abu_filt <- st_len_filt[,abundances >= MINABUND]


## ----compare filter steps 2 -- not working------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


## ----tibbles tracking reads through pipeline----------------------------------

track2 <- as_tibble(track)
glimpse(rownames(track2))

track2 %>% 
  summarize_all(median)
  
# Write track
write.csv(track, file = "track.csv")



## ----seqs to FASTA------------------------------------------------------------
getSequences(seqtab.nochim) 

#export seqs as fasta
uniquesToFasta(getUniques(seqtab.nochim), fout="./uniqueSeqs_DINO.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab.nochim)))))

## confirmed primers completed removed! 

## Original code: uniquesToFasta(getUniques(seqtab.nochim), fout="~/Dropbox/1_symITSdata/uniqueSeqs.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab.nochim))))) 


## ----decipher phylogenetic tree-----------------------------------------------
# alignment using Decipher package note: may need to load
# This is for the length and abundance filtered samples.
seqs_filt <- getSequences(st_len_abu_filt)
names(seqs_filt) <- seqs_filt # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs_filt), anchor=NA)

## AlignSeqs : Align A Set Of Unaligned Sequences - Performs profile-to-profile alignment of multiple unaligned sequences following a guide tree.



## -----------------------------------------------------------------------------
uniquesToFasta(getUniques(st_len_abu_filt),fout="./uniqueSeqs.fasta")

#export seqs as fasta
uniquesToFasta(getUniques(st_len_abu_filt),fout= "./uniqueSeqs.fasta", ids=paste0("asv", seq(length(getUniques(st_len_abu_filt)))))


## -----------------------------------------------------------------------------
#used phanforn to make tree

## Change sequence alignment output into a phyDat structure 

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

## distance matrix
dm <- dist.ml(phang.align)

## Unweighted pair group method with arithmetic mean

treeUPGMA <- upgma(dm)

### In bioinformatics, UPGMA is used for the creation of phenetic trees (phenograms). UPGMA was initially designed for use in protein electrophoresis studies, but is currently most often used to produce guide trees for more sophisticated algorithms. This algorithm is for example used in sequence alignment procedures, as it proposes one order in which the sequences will be aligned. Indeed, the guide tree aims at grouping the most similar sequences, regardless of their evolutionary rate or phylogenetic affinities, and that is exactly the goal of UPGMA

### In phylogenetics, UPGMA assumes a constant rate of evolution (molecular clock hypothesis) and that all sequences were sampled at the same time, and is not a well-regarded method for inferring relationships unless this assumption has been tested and justified for the data set being used. Notice that even under a 'strict clock', sequences sampled at different times should not lead to an ultrametric tree

## The most important practical issues: UPGMA provides rooted tree as a result, while NJ unrooted, and you have to take care proper rooting the NJ tree afterward. Also, UPGMA is regarded as an unreliable method, so I would prefer to use NJ (see ultrametric discussion above)

## UPGMA: equal rate among branches >> equal depth from root ; NJ: allow diff rate among branches >> diff depth from root

## Neighbor joining 
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## plot nj tree

layout(matrix(c(1,2), 2, 1), height=c(1,2))
par(mar = c(0,0,2,0)+ 0.1)

# plot(treeUPGMA, main="UPGMA")

plot(treeNJ, "unrooted", main="NJ") 
fit = pml(treeNJ, data=phang.align)

# NJ pref of tree -- these trees go into phyloseq (collapse, plot based on genetic distance)

write.tree(treeNJ, file="nj.tre")

# write.tree(treeUPGMA, file="upgma.tre")

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

detach("package:phangorn", unload=TRUE)


## -----------------------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())


## ----handoff to Phyloseq------------------------------------------------------

#Import sample_data

samdf <- read.csv("/home/vmglynn/projects/def-barrett/vmglynn/internship2020/FULL_DINO/CSVver5_Poc_ITS_SeqInfo_3May2021_EnvData.csv", row.names = 1)

stab <- st_len_abu_filt

# Code to make table based on sample name, if informative 



## ----Phyloseq object----------------------------------------------------------
ps <- merge_phyloseq(otu_table(stab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               #tax_table(taxa_filt),
               phy_tree(treeNJ))

## sample_data(samdf2) -- sample_data =  st_len_abu_filt ?


ps
saveRDS(samdf, "samdf.rds")
write.csv(samdf, file = "samdf.csv")
#saveRDS(taxa_filt, "taxa_filt.rds")
saveRDS(treeNJ, "treeNJ.rds")
saveRDS(stab, "stab.rds")
saveRDS(ps, "symITSps.rds")

## ----phyloseq ExportFASTA-----------------------------------------------------


# export_fasta(ps = symITSps.rds, file = NULL, rank = NULL)

readRDS("symITSps.rds")


