#!/usr/bin/env Rscript

# ------------------------------------------------------------------------
# Preparation
# ------------------------------------------------------------------------

#clear the object from memory
rm(list=ls())

#### CHANGE THIS PATHS ####
source <- "/scicore/home/schlae0003/GROUP/projects/P2_Jan_microbials_experiment/roots/02_ASV"
path.in <- "../01_demultiplex/out/" #needs a "/" at the end
path.primer <- "../01_demultiplex/" #needs a "/" at the end
path.rds <- "RDS/"
path.out <- "output/"
taxa_database <- "../../../../taxanomy_databases/utax_reference_dataset_10.05.2021.fasta"

##########################

setwd(source)

#libaries
library(dada2)
library(Biostrings)
library(ShortRead)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(phyloseq)
library(xlsx)
library(seqinr)
library(vegan)
library(tibble)
library(dplyr)
library(DECIPHER)
library(Biostrings)
library(qdapRegex)

#generate output files
output <- T

#list of bacteria runs
runs <- list.files(path.in)

#create folders
dir.create(path.rds)
dir.create(path.out)

#settings
rc <- dada2:::rc
theme_set(theme_bw())

# ------------------------------------------------------------------------
# Loop over all bacteria runs
# ------------------------------------------------------------------------
for (run in runs) {
  print(paste(run, "started"))
  
  # ------------------------------------------------------------------------
  # DADA2 pipeline - PacBio
  # ------------------------------------------------------------------------
  
  #file_names
  path_run <- paste(path.in, run, sep="")
  fns <- list.files(path_run, pattern="fastq", full.names=TRUE)
  assign(paste(run, "_fns", sep=""), fns)
  
  #get primers
  path_primers <- paste(path.primer, run, "/primer.fasta", sep="")
  primers <- read.fasta(path_primers, as.string = T, forceDNAtolower = F)
  primerF <- primers[["primerF"]][1]
  primerR <- primers[["primerR"]][1]

  #remove primers
  nop <- paste0(path.in, run, "/noprimers/", basename(fns))
  prim <- removePrimers(fns, nop, primer.fwd = primerF, primer.rev = dada2::rc(primerR) , orient =T, verbose =T)
  
  #sequence length before filtering
  lens.fn <- lapply(fns, function(fn) nchar(getSequences(fn)))
  lens <- do.call(c, lens.fn)
  histo_unfiltered <- hist(lens, 100, main = paste(run, "unfiltered"))
  assign(paste(run, "histo_unfiltered", sep="_"), histo_unfiltered)
  
  #names filter
  filts <- file.path(path_run, "filtered", basename(fns))
  assign(paste(run, "_filts", sep=""), filts)

  #filter
  #max Len prevents double and triple CCS reads to pass the filter
  track <- filterAndTrim(nop, filts, minQ=3, minLen=500, maxLen=1800, maxN=0, rm.phix=T, maxEE=2, multithread = T)

  assign(paste(run, "track", sep="_"), track)
  saveRDS(track, paste(path.rds, run, "_out.RDS", sep=""))
  
  #sequence length after filtering
  histo_filtered <- hist(lens, 100, main = paste(run, "filtered"))
  assign(paste(run, "histo_filtered", sep="_"), histo_filtered)
  
  #dereplicate
  drp <- derepFastq(filts, verbose=TRUE)
  assign(paste(run, "drp", sep="_"), drp)
  saveRDS(drp, paste(path.rds, run, "_derep.RDS", sep=""))
  
  #error estimation for PacBio
  err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE, qualityType="FastqQuality")
  assign(paste(run, "err", sep="_"), err)
  saveRDS(err, paste(path.rds, run, "_err.RDS", sep=""))
  
  #denoise with inference algorithm
  dd <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE, selfConsist=TRUE)
  assign(paste(run, "dd", sep="_"), dd)
  saveRDS(dd, paste(path.rds, run, "_denoise.RDS", sep=""))
  
}

# ------------------------------------------------------------------------
#combine all runs
# ------------------------------------------------------------------------

#combine and save
fns <- vector()
for(i in runs){fns <- c(fns, get(paste(i, "_fns", sep="")))}
saveRDS(fns, paste(path.rds, "fns.RDS", sep=""))

filts <- vector()
for(i in runs){filts <- c(filts, get(paste(i, "_filts", sep="")))}
saveRDS(filts, paste(path.rds, "filts.RDS", sep=""))

track <- get(paste(runs[1], "_track", sep=""))
if(length(runs)>1){for(i in runs[2:length(runs)]){track <- rbind(track, get(paste(i, "_track", sep="")))}}
saveRDS(track, paste(path.rds, "track.RDS", sep=""))

drp <- vector()
for(i in runs){drp <- c(drp, get(paste(i, "_drp", sep="")))}
saveRDS(drp, paste(path.rds, "drp.RDS", sep=""))

dd <- vector()
for(i in runs){dd <- c(dd, get(paste(i, "_dd", sep="")))}
saveRDS(dd, paste(path.rds, "dd.RDS", sep=""))

#names from all bacterial runs
sample.names <- sapply(strsplit(basename(fns), "_"), `[`, 2)

#sequence table
seqtab <- makeSequenceTable(dd)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#dim(seqtab.nochim)

#tracking
getN <- function(x) sum(getUniques(x))
track <- cbind(track, sapply(dd, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "nochim")
saveRDS(track, paste(path.rds, "track.RDS", sep=""))

#taxa
taxa <- assignTaxonomy(seqtab.nochim, taxa_database, multithread=TRUE,tryRC = TRUE ) #try in both directions

#taxa confidence intervall
# taxa_conf <- assignTaxonomy(seqtab.nochim, taxa_database, multithread=TRUE,tryRC = TRUE, minBoot = 0, outputBootstraps = TRUE) #tries in both directions

#change headers from the long sequence to simple ASV numbers
seq <- colnames(seqtab.nochim)
ASV <- c()
for (i in 1:length(seq)) {ASV[i] <- paste("ASV",i , sep="") }
#ASV_seq <- cbind(ASV, seq)
#colnames(ASV_seq) <- c("ASV", "seq")
ASV_seq <- data.frame(ASV=ASV, seq=seq)
ASV_seq[,1] <- paste(">", ASV_seq[,1], sep="")
saveRDS(ASV_seq, paste(path.rds, "ASV_seq.RDS", sep=""))

#replace sequence headers by ASV numbers in seqtab.nochim
colnames(seqtab.nochim) <- ASV
saveRDS(seqtab.nochim, paste(path.rds, "seqtab.nochim.RDS", sep=""))

#replace the sequence headers by the correct ASV number in taxa
for (i in 1:nrow(taxa)){
  for (j in 1:nrow(ASV_seq)){
    if(rownames(taxa)[i] == ASV_seq[j,2]){
      rownames(taxa)[i] <- ASV_seq[j,1]
    }
  }
}
rownames(taxa) <- gsub(">", "", rownames(taxa))

#"UTAX db" has different FASTA headers than "general release db"
#taxomony df needs to be done manually

#empty df
temp <- data.frame(matrix(NA, nrow(taxa), 7))
rownames(temp) <- rownames(taxa)
colnames(temp) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species") 

#extract taxa
for(tax in 1:nrow(taxa)){
  tax_string <- paste0(taxa[tax,2],",") #extract string and add "," at the end
  tax_vector <- qdapRegex::ex_between(tax_string, ":", ",")[[1]] #use regex to extract taxa
  temp[tax, 1:length(tax_vector)] <- tax_vector #add taxa
}

#remove db-number
temp[,7] <- gsub("_SH[0-9]*.08FU", "", temp[,7])

#overwrite taxa df
taxa <- temp; rm(temp)
saveRDS(taxa, paste(path.rds, "taxa.RDS", sep=""))
  
# ------------------------------------------------------------------------
# Clustering
# ------------------------------------------------------------------------

nproc <- 16

#change ASV names to corresponding sequences
seqtab.nochim_ASV <- seqtab.nochim
colnames(seqtab.nochim_ASV) <- ASV_seq$seq

#align sequences and compute distance matrix
ASV_seqs <- Biostrings::DNAStringSet(ASV_seq$seq) #convert to DNA string set
ASV_seqs_aln <- DECIPHER::AlignSeqs(ASV_seqs, processors = nproc)
ASV_dist <- DECIPHER::DistanceMatrix(ASV_seqs_aln, processors = nproc)

#OTUs
ASV97_cluster <- DECIPHER::IdClusters(ASV_dist, method = "complete", processors = nproc, cutoff = 0.03) #cutoff = 0.03 -> 97% OTU 
ASV98_cluster <- DECIPHER::IdClusters(ASV_dist, method = "complete", processors = nproc, cutoff = 0.02) #cutoff = 0.02 -> 98% OTU 
ASV99_cluster <- DECIPHER::IdClusters(ASV_dist, method = "complete", processors = nproc, cutoff = 0.01) #cutoff = 0.01 -> 99% OTU

#loop over the 3 thresholds
for (n in c("97", "98", "99")) {

    #get the data
  ASV_cluster <- get(paste0("ASV", n, "_cluster"))
  ASV_cluster$ASV <- gsub(">","",ASV_seq$ASV) #add ASV name
  ASV_cluster$ASV_seq <- ASV_seq$seq #add sequences
  
  # Number of sequences per ASV
  ASV_cluster$ASV_abu <- NA
  for (seq in ASV_cluster$ASV_seq) {
    abu <- sum(seqtab.nochim_ASV[,colnames(seqtab.nochim_ASV) == seq])
    ASV_cluster$ASV_abu[ASV_cluster$ASV_seq == seq] <- abu
  }
  
  #most abundant sequences for each OTU
  ASV_cluster$OTU_seq <- NA
  ASV_cluster$OTU_abu <- NA
  
  for (OTU in ASV_cluster$cluster) {
    ASV_cluster_subset <- ASV_cluster[ASV_cluster$cluster == OTU,] #all entries for current OTU
    n_seq <- sum(ASV_cluster_subset$ASV_abu)
    most_abu_seq <- ASV_cluster_subset$ASV_seq[ASV_cluster_subset$ASV_abu == max(ASV_cluster_subset$ASV_abu)]# OTU(s) with most seqs
    most_abu_seq <- most_abu_seq[1] #take the first entry in case multiple OTUs have max number of seqs
    ASV_cluster$OTU_seq[ASV_cluster$cluster == OTU] <- most_abu_seq #add seq to table
    ASV_cluster$OTU_abu[ASV_cluster$cluster == OTU] <- n_seq #add num to table
  }
  
  #sort by OTU abundance
  ASV_cluster_sort <- ASV_cluster[order(ASV_cluster$OTU_abu, decreasing = T),] #sort by OTU name
  
  #rename OTU by abundance
  clusters_sort <- unique(ASV_cluster_sort$cluster)
  ASV_cluster$OTU <- NA
  for(i in seq(length(clusters_sort))){
    ASV_cluster$OTU[ASV_cluster$cluster == clusters_sort[i]] <- paste0("OTU", i) #rename
  }
  
  ### FINAL OTU FILES
  
  #OTU ASV_seq
  OTU_ASV_seq_all <- ASV_cluster[order(ASV_cluster$OTU_abu, decreasing = T),]
  OTU_ASV_seq_all <- data.frame(OTU = OTU_ASV_seq_all$OTU, seq = OTU_ASV_seq_all$OTU_seq)
  OTU_ASV_seq <- data.frame(OTU = unique(OTU_ASV_seq_all$OTU), seq = NA)
  for (OTU in OTU_ASV_seq$OTU) {
    seq <- OTU_ASV_seq_all$seq[OTU_ASV_seq_all$OTU == OTU][1] #get the sequence
    OTU_ASV_seq$seq[OTU_ASV_seq$OTU == OTU] <- seq #add to seqtab.nochima frame
  }
  
  
  #OTU seqtab.nochim
  OTU_seqtab.nochim <- seqtab.nochim_ASV %>% t %>% rowsum(ASV_cluster$cluster) %>% t #cluster
  for(i in seq(length(clusters_sort))){
    colnames(OTU_seqtab.nochim)[colnames(OTU_seqtab.nochim) == clusters_sort[i]] <- paste0("OTU", i) #rename
  }
  OTU_seqtab.nochim <- OTU_seqtab.nochim[,OTU_ASV_seq$OTU] #sort
  
  
  #OTU taxa
  OTU_taxa <- data.frame(matrix(nrow = nrow(OTU_ASV_seq), ncol = 7))
  rownames(OTU_taxa) <- OTU_ASV_seq$OTU
  colnames(OTU_taxa) <- colnames(taxa)
  
  for(OTU in OTU_ASV_seq$OTU){
    seq <- OTU_ASV_seq$seq[OTU_ASV_seq$OTU == OTU] #OTU seq
    ASV <- ASV_cluster$ASV[ASV_cluster$ASV_seq == seq] #ASV belong to OTU seq

    taxa_ASV <- taxa[rownames(taxa) == ASV[1], ] #taxa belong to the first ASV
    OTU_taxa[rownames(OTU_taxa) == OTU, ] <- taxa_ASV #transfer to OTU taxa
  }
  
  #OTU_ASV
  OTU_ASV <- data.frame(OTU = ASV_cluster$OTU, ASV = ASV_cluster$ASV, 
                        OTU_seq = ASV_cluster$OTU_seq, ASV_seq = ASV_cluster$ASV_seq,
                        OTU_abu = ASV_cluster$OTU_abu, ASV_abu = ASV_cluster$ASV_abu)
  
  OTU_ASV <- OTU_ASV[order(ASV_cluster$OTU_abu, decreasing = T),] #sort
  
  #store clustered data
  assign(paste0("OTU_ASV_seq", n), OTU_ASV_seq)
  assign(paste0("OTU_seqtab.nochim", n), OTU_seqtab.nochim)
  assign(paste0("OTU_taxa", n), OTU_taxa)
  assign(paste0("OTU_ASV", n), OTU_ASV)
  
  #save RDS files
  dir.create(paste0(path.rds,"ASV",n))
  saveRDS(OTU_ASV_seq, paste0(path.rds,"ASV",n, "/SEQ", n, ".RDS"))
  saveRDS(OTU_seqtab.nochim, paste0(path.rds,"ASV",n, "/DAT", n, ".RDS"))
  saveRDS(OTU_taxa, paste0(path.rds,"ASV",n, "/TAXA", n, ".RDS"))
  saveRDS(OTU_ASV, paste0(path.rds,"ASV",n, "/ABU", n, ".RDS"))
}
  
# ------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------
if(output){
  
  #QUALITY CHECK
  
  dir.create(paste0(path.out,"quality_check"))
  
  #number of reads
  write.xlsx(track, paste0(path.out, "quality_check/track.xlsx"))
  
  #length histo
  pdf(paste0(path.out,"quality_check/reads_length.pdf"))
    for(i in runs){
      plot(get(paste(i, "histo_unfiltered", sep="_")), main = paste(i, "unfiltered"))
      plot(get(paste(i, "histo_filtered", sep="_")), main = paste(i, "filtered"))
    }
  dev.off()
  
  #quality profiles unfiltered and untrimmed
  pdf(paste0(path.out,"quality_check/reads_quality_unfilt_untrim.pdf"))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(fns[i]))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()
  
  #quality profiles filtered and trimmed
  pdf(paste0(path.out,"quality_check/reads_quality_filt_trim.pdf"))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(filts[i]))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()
  
  #error rate learning
  runs_err <- c (paste0(runs, "_err"))
  
  pdf(paste0(path.out,"quality_check/error_rate.pdf"))
  for (i in 1:length(runs)) {
    plot(plotErrors(get(runs_err[i]), nominalQ=TRUE))
    grid.text(runs[i], hjust=-1.5, vjust = -27.5, rot = 90)
  }
  dev.off() 
  
  #rarefaction
  pdf(paste0(path.out,"quality_check/rarefaction.pdf"))
  for (run in runs) {
    #sp <- paste0(substr(run,1,1), substr(run, nchar(run), nchar(run))) #search pattern
    sp <- gsub("microbials", "m", run) #search pattern
    DAT <- seqtab.nochim[grep(sp, rownames(seqtab.nochim)),]
    rarecurve(DAT, step = 20, label=F, main = paste("Rarefaction", run), ylab="ASVs", xlab="Sequencing depth")
  }
  dev.off()
  
  print("quality check output done")
  
  # ASV TABLES
  
  #97
  dir.create(paste0(path.out,"ASV97"))
  write.table(OTU_seqtab.nochim97, paste0(path.out, "ASV97/DAT97.tab"), sep="\t")
  write.table(OTU_taxa97, paste0(path.out, "ASV97/TAXA97.tab"), sep="\t")
  write.table(OTU_ASV_seq97, paste0(path.out, "ASV97/SEQ97.tab"), sep="\t", row.names = F)
  write.table(OTU_ASV97, paste0(path.out, "ASV97/ABU.tab"), sep="\t", row.names = F)
  
  #98
  dir.create(paste0(path.out,"ASV98"))
  write.table(OTU_seqtab.nochim98, paste0(path.out, "ASV98/DAT98.tab"), sep="\t")
  write.table(OTU_taxa98, paste0(path.out, "ASV98/TAXA98.tab"), sep="\t")
  write.table(OTU_ASV_seq98, paste0(path.out, "ASV98/SEQ98.tab"), sep="\t", row.names = F)
  write.table(OTU_ASV98, paste0(path.out, "ASV98/ABU.tab"), sep="\t", row.names = F)
  
  #99
  dir.create(paste0(path.out,"ASV99"))
  write.table(OTU_seqtab.nochim99, paste0(path.out, "ASV99/DAT99.tab"), sep="\t")
  write.table(OTU_taxa99, paste0(path.out, "ASV99/TAXA99.tab"), sep="\t")
  write.table(OTU_ASV_seq99, paste0(path.out, "ASV99/SEQ99.tab"), sep="\t", row.names = F)
  write.table(OTU_ASV99, paste0(path.out, "ASV99/ABU.tab"), sep="\t", row.names = F)
  
  #100
  dir.create(paste0(path.out,"ASV100"))
  write.table(seqtab.nochim, paste0(path.out, "ASV100/DAT100.tab"), sep="\t")
  write.table(taxa, paste0(path.out, "ASV100/TAXA100.tab"), sep="\t")
  #write.table(taxa_conf, paste0(path.out, "ASV100/TAXA100_conf.tab"), sep="\t")
  write.table(ASV_seq, paste0(path.out, "ASV100/SEQ100.tab"), sep="\t")
  
}
