library(dada2)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(dplyr)
library(Biostrings)
library(gridExtra)
library(DECIPHER)
library(xlsx)
library(ape)
theme_set(theme_bw())

path <- "C:\\~directory containing fastq files"
pathtab <- paste0("C:\\~directory to save .csv files")
list.files(path)
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])  
plotQualityProfile(fnRs[1:2])  
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)  
head(out)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
file.exists(fnFs)
file.exists(fnRs)
errF <- learnErrors(filtFs, multithread=TRUE)  
errR <- learnErrors(filtRs, multithread=TRUE)  
plotErrors(errF, nominalQ=TRUE)  
derepFs <- derepFastq(filtFs, verbose=TRUE) 
derepRs <- derepFastq(filtRs, verbose=TRUE) 
names(derepFs) <- sample.names # Name the derep-class objects by the sample names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE) 
dadaRs <- dada(derepRs, err=errR, multithread=TRUE) 
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
length(mergers)
head(mergers[[1]])
seqtabe <- makeSequenceTable(mergers)
dim(seqtabe)
table(nchar(getSequences(seqtabe)))
plot(table(nchar(getSequences(seqtabe))))
seqtab <- seqtabe[,nchar(colnames(seqtabe)) %in% seq(300,405)]
dim(seqtab)
table(nchar(getSequences(seqtabe)))
plot(table(nchar(getSequences(seqtab))), xlab = "Reads R1+R2 merged length")
#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtabe, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtabe)

sequences <- getSequences(seqtab.nochim)
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
#tracking reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


#assigning taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr99_v138_train_set.fa.gz", multithread=FALSE)
taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

setwd(pathtab)
file.remove(list.files())
setwd(path)
write.csv(t(seqtab.nochim), paste0(pathtab, "seqtab-nochim.csv"), quote=FALSE)
seqs <- colnames(seqtab.nochim)
SeqName <- vector(dim(seqtab.nochim)[2], mode="character")
SeqName_ft <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  SeqName[i] <- paste("seq", i, sep="")
  SeqName_ft[i] <- paste(">seq", i, sep="")
}
fastaseqs_ft <- rbind(SeqName_ft, seqs)
write(fastaseqs_ft, paste0(pathtab,"fastaseqs.fasta"))
fastaseqs <- cbind(SeqName, seqs)
write.csv(fastaseqs, paste0(pathtab,"fastaseqs.csv"), quote=FALSE, row.names = FALSE)
seqtab.nochim.t <- t(seqtab.nochim)
row.names(seqtab.nochim.t) <- SeqName
seqtab.nochim.t <- tibble::rownames_to_column(as.data.frame(seqtab.nochim.t), "SeqName")
write.csv(seqtab.nochim.t, paste0(pathtab, "counts.csv"), quote=FALSE, row.names = FALSE)
taxtable <- taxa
row.names(taxtable) <- SeqName
taxtable <- tibble::rownames_to_column(as.data.frame(taxtable), "SeqName")
write.csv(taxtable, paste0(pathtab, "taxtable.csv"), quote=FALSE, row.names = FALSE)
tax_counts <- left_join(taxtable, seqtab.nochim.t, by = "SeqName", keep = TRUE)
tax_counts_fasta <- left_join(tax_counts, as.data.frame(fastaseqs), by = c("SeqName.x" =
                                                                             "SeqName"), keep = TRUE)
openxlsx::write.xlsx(tax_counts_fasta, file = paste0(pathtab, "tax_counts_fasta.xlsx"),
                     overwrite = TRUE, asTable = FALSE, sheetName = "EAW_18S", firstRow = TRUE, zoom = 90,
                     keepNA = TRUE)
save.image(paste0("C:\\~full path of location"))

path <- ("C:\\~pathtab")
#seq_table <- read.csv("otumat.csv")
taxa_matrix <- read.csv("tablescounts.csv")
otumat <- as.matrix(taxa_matrix)
print(otumat)
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
otumat
#sampledata1 <- read.table("metadata1.txt", header = TRUE, row.names=1, fill=TRUE, sep = "\t")
taxa_table <- read.csv("tablestaxtable.csv")
taxmat <- as.matrix(taxa_table)
print(taxmat)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
taxmat
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
TAX
physeq = phyloseq(OTU, TAX)
physeq
plot_bar(physeq, fill = "Family")
plot_bar(physeq, fill = "Genus")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)
physeq1 = merge_phyloseq(physeq, random_tree)
physeq1
sampledata <- physeq1
physeq2 = phyloseq(OTU, TAX, sampledata, random_tree)
physeq2
identical(physeq1, physeq2)
plot_tree(physeq1, color="Location", label.tips="taxa_names", ladderize="left", plot.margin=0.3)