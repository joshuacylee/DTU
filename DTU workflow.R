#Before begining this pipeline, do make sure working directory contains the following:
#Salmon folder, with subfolder containing NAME of the sample, and quant.sf in the subfolder
#Need to create samples.txt with details of samples that are being compared. Format explained below

#This script uses recurrent vs non-recurrent ccRCC samples as an example

#load libraries
library(GenomicFeatures)
library(DRIMSeq)
library(tximport)
library(rnaseqDTU)
library(DESeq2)
library(stageR)
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(readr)
library(biomaRt)

#Set up sample table (from txt file) which consists of 2 columns:
#First column: sample_id. Try not to start with number and sample_id MUST match folder name containing salmon quant file
#Second column: condition. It is possible to have multiple comparisons, but here only cancer recurrence status are used.

samps <- read.table("samples.txt", header = TRUE, sep = '\t')

#check if columns names, and if contents look ok
head(samps)

#set up factors (variable for comparisons - important for much later in the pipeline!), and optionally order by wt vs ko, or naive vs th1 etc. Match with what's in samples.txt
samps$condition <- factor(samps$condition)
samps$condition <- ordered(samps$condition, levels = c("nonrecurrent", "recurrent"))
table(samps$sample_id)
table(samps$condition)

#set up file paths - 1 level up from the folder names appeared in samples.txt
files <- file.path("/file/to/path/salmon_quant", samps$sample_id, "quant.sf")
names(files) <- samps$sample_id
head(files)

#check if all files exist. If TRUE, proceed.
all(file.exists(files))

#make a table with all ensembl genes vs corresponding transcript ID. This allows us to assign gene IDs to each transcripts later on

txdb <- makeTxDbFromBiomart(biomart="ensembl",dataset = "hsapiens_gene_ensembl")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)


#To avoid differences in ID version number (the .X after the ENST/ENSG), delete them from tx2gene object. change '18' accordingly - essentially to get rid of the isoform version
##This is optional since tximport has an ignoreTxVersion step following on
tx2gene$TXNAME <- substr(tx2gene$TXNAME, 0, 15)
tx2gene$GENEID <- substr(tx2gene$GENEID, 0, 15)
head(tx2gene)

#import salmon quant files. countsFromAbundance="dtuscaledTPM" <- In the publication they recommend scaledTPM but since then the authors have introduced dtuscaledTPM without updating the paper.
#ignoreAfterBar = TRUE is to avoid taking in all info from genecode gtf file. ignoreTxVersion = TRUE here - some versions this does not work. If so make sure there's no isoform version following on, using previous step.
txi.salmon <- tximport(files, type="salmon", txOut=TRUE, tx2gene = tx2gene,
                       countsFromAbundance="dtuscaledTPM", ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)
head(txi.salmon)

#get rid of all 0 counts transcripts (ones that are not expressed at all across samples, hence rowSums) and produce a csv file with dtuscaledTPM for plotting
##Requirement for number of reads can be adjusted later on in the pipeline.
cts <- txi.salmon$counts
cts <- cts[rowSums(cts) > 0,]
write.csv(cts, file="cts.csv")

#read file and add columns, then delete the version number. Write over the previously saved file.
file = read.csv(file = "cts.csv", header = TRUE, sep = ",")
head(file)
colnames(file)[1] <- "TXNAME"
file$TXNAME <- substr(file$TXNAME, 0, 15)
head(file)

#making sure all mapped items from file are present in tx2gene - make note of number/% of dropped out transcripts.
#This is due to the differences between reference transcriptome used.
file <- file[file$TXNAME %in% tx2gene$TXNAME, ]
all(file$TXNAME %in% tx2gene$TXNAME)
cts <- file
write.csv(file, file="cts.csv")

#if index column appears, can use cts <- cts[,-1] to delete column

#Add corresponding geneID and symbols to table using biomaRT package. useEnsembl dataset can be changed for different species (hsapiens_gene_ensembl)
#Ensembl site can be unresponsive sometimes - try again if it is not working, just to make sure the objects exist!
mart = useMart('ensembl')
ensembl = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

##In case connection fails multiple times, chances are getBM won't work in this mirror. Try this instead (uswest or asia or useast):
## Mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirroe = "uswest")
### And instead of mart = ensembl below, use mart = mart.

#add gene id and gene name to the list. Can add more attributes here if they are needed. If txname geneID requires version, add "_version"after "id". values = file$TXNAME)
gene <- getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id","external_gene_name"),values = file$TXNAME, mart = ensembl)
head(gene)

#match lists and see if all values exist in your file
id <- match(file$TXNAME , gene$ensembl_transcript_id)

#add matched columns (geneID and symbol for transcripts > 0 in the samples)
file$gene_id <- gene$ensembl_gene_id[id]
head(file)

#Do it again for gene symbol
id2 <- match(file$gene_id, gene$ensembl_gene_id)
file$gene_symbol <- gene$external_gene_name[id2]
head(file)

#move columns to make it them look more organised
file = file %>% relocate(gene_id, .after = TXNAME)
file = file %>% relocate(gene_symbol, .after = gene_id)

#export file as final table for scaledTPM
write.csv(file, file = "ccRCC_dtuscaledTPM.csv")

#Going back to previous files for DRIMseq input, now we have to build a txdf object (tx to g dataframe)
#make sure number of txdf rows match up with cts
txdf <- tx2gene[match(cts$TXNAME,tx2gene$TXNAME),]
head(txdf)
all(cts$TXNAME %in% txdf$TXNAME)

#If TRUE, build count table
counts <- data.frame(gene_id=txdf$GENEID,
                     feature_id=txdf$TXNAME,
                     cts)

head(counts)

#start DRIMseq analysis - make note of how many genes are detected at this stage. Keep in mine these are none-0 counts, i.e. mapped in at least 1 sample
d <- dmDSdata(counts=counts, samples=samps)
d

#pulling out first row to see what can be seen
methods(class=class(d))
counts(d[1,])[,1:8]

#check number of transcripts mapped per gene
table(table(counts(d)$gene_id))

#define sample groups, n being all the samples and n.small being number of samples within condition
#Filter low abundant transcripts- total number of samples = n = 12, each group = n.small = 6
#rule: we want at least a count of 10 mapped transcripts in at least number of samples within small group (in this case:5)
#also the relative abundance proportion of at least 0.05 (5%) in at least n.small samples
#ALL samples must express the gene (not necessarily the specific transcript), with at least 5 here.
##This can be adjusted depending on stringency
n <- 12
n.small <- 6
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=5,
              min_samps_feature_prop=n.small, min_feature_prop=0.05,
              min_samps_gene_expr=n, min_gene_expr=5)


#check number of genes still remain
d
table(table(counts(d)$gene_id))

#design experiment - keep note of what colnames are - first colname:Intercept, second colname is crucial for next step
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)

#check dispersion (long wait). For coef = "", where ""is = the second factor of colnames(design_full) <- PLEASE MAKE SURE THEY MATCH
#Save data at this point to avoid crash and lost of data
set.seed(1)
system.time({
  d <- dmPrecision(d, design=design_full)
  d <- dmFit(d, design=design_full)
  d <- dmTest(d, coef="")
})

#produce DRIMseq results table for further quantification
res <- DRIMSeq::results(d)
head(res)

#Add transcripts ID (ENST) back in
res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)

#Eliminate all NAs from p value
no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)

#add gene symbols back into the table
id3 <- match(res.txp$gene_id, gene$ensembl_gene_id)
res.txp$gene_symbol <- gene$external_gene_name[id3]
res.txp = res.txp %>% relocate(gene_symbol, .after = feature_id)


#export res.txp file -> DRIMseq results with pvalue. Note that this is before the 2nd stage statistical filtering by stageR,
#Useful table to plot proportions ensembl_transcript_id
write.csv(res.txp, file="ccRCC_DRIMseq_padj.csv")
idx <- which(res$adj_pvalue < 0.1)[1]
res[idx,]

#StageR Stats -> essentially takes it a step further:
#DRIMSeq looks at whether proportion changes are significant, where as StageR look at whether transcripts or gene changes significantly, and must satisfy both to be significant

nrow(res)
nrow(res.txp)
pScreen <- res$pvalue
names(pScreen) <- sub("(\\.)","p",res$gene_id,perl=TRUE)
pConfirmation <- matrix(res.txp$pvalue, ncol=1)
rownames(pConfirmation) <- gsub("(\\.)","p",res.txp$feature_id,perl=TRUE)
tx2gene <- res.txp[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- gsub("(\\.)","p",tx2gene[,i],perl=TRUE)



#At this point we can decide on what alpha/pvalue denotes stats sig. Here, 0.1 is used.
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=TRUE,
                                  onlySignificantGenes=TRUE)
})
head(drim.padj)

#add gene symbols back into the table
id4 <- match(drim.padj$geneID, gene$ensembl_gene_id)
drim.padj$gene_symbol <- gene$external_gene_name[id4]
drim.padj = drim.padj %>% relocate(gene_symbol, .after = geneID)

#export significance table
write.csv(drim.padj, file="ccRCC_DRIMseq_StageR_padj.csv")

#convert counts table back to usable format for the next few steps and caluclate isoform proportions
#export proportions - only take genes which are listed from drim.padj

counts=dplyr::select(counts,-TXNAME)
drim.prop = reshape2::melt(counts[counts$feature_id %in% proportions(d)$feature_id,], id = c("gene_id", "feature_id"))
drim.prop = drim.prop[order(drim.prop$gene_id, drim.prop$variable, drim.prop$feature_id),]

system.time({
  drim.prop = drim.prop %>%
    group_by(gene_id, variable) %>%
    mutate(total = sum(value)) %>%
    group_by(variable, add=TRUE) %>%
    mutate(prop = value/total)
})

drim.prop = reshape2::dcast(drim.prop[,c(1,2,3,6)], gene_id + feature_id ~ variable)
id5 <- match(drim.prop$gene_id, gene$ensembl_gene_id)
drim.prop$gene_symbol <- gene$external_gene_name[id5]
drim.prop = drim.prop %>% relocate(gene_symbol, .after = feature_id)
write.csv(drim.prop, file = "ccRCC_DRIMSeq_prop.csv")

#DRIMseq proportions of genes that are marked as significant from drim.padj
drim.prop.sig = drim.prop[drim.prop$gene_id %in% drim.padj$geneID,]
id6 <- match(drim.padj$geneID, drim.prop.sig$gene_id)
drim.prop.sig$gene_padj <- drim.padj$gene[id6]
drim.prop.sig = drim.prop.sig %>% relocate(gene_padj, .after = gene_symbol)
id7 <- match(drim.padj$txID, drim.prop.sig$feature_id)
drim.prop.sig$transcript_padj <- drim.padj$transcript[id7]
drim.prop.sig = drim.prop.sig %>% relocate(transcript_padj, .after = gene_padj)
write.csv(drim.prop.sig, file = "ccRCC_DRIMSeq_prop_sig.csv")

#at this point we can plot transcripts proportion in R (although with the csv file it's probably better to plot of graphpad)
gene_id = "ENSG" #substitute ENSG gene ID for gene of interest
plotProportions(d, gene_id = gene_id, group_variable = "condition")

#DEXseq also utilises data processed from the dmPrecision step, and after that filtered by StageR
##DEXSeqDataSet doesn't like having characters at the header for samps - DRIMSeq::samples(d) for sample.data doesn't work. Input from file again with row.names = 1

sample.data <- read.table("samples.txt", header = TRUE, sep = '\t', row.names = 1)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)

#check dxd file
dxd

#calculate size factors and dispersions
system.time({
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)
})

#Calculate per gene/ per transcript adjusted p value -
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval = perGeneQValue(dxr)
dxr.g = data.frame(gene = names(qval), qval)
dxr.t = as.data.frame(dxr[, c("featureID","groupID","pvalue")])
head(dxr.t)
write.csv(dxr.g, file="ccRCC_DEXseq_genes.csv")
write.csv(dxr.t, file="ccRCC_DEXseq_transcripts.csv")

#Start stageR - strp function strips dictates 15 letters ENST/ENSG to make sure there's no version number - change if length is different (i.e human)
strp <- function(x) substr(x,1,15)
pScreen = qval
names(pScreen) = strp(names(pScreen))
pConfirmation = matrix(dxr.t$pvalue, ncol=1)
dimnames(pConfirmation) = list(strp(dxr.t$featureID),"transcript")

#set of genes/trasncripts passing screening with padj less than 510% - change alpha(padj) if need be
stageRObj = stageRTx(pScreen = pScreen,
                     pConfirmation = pConfirmation,
                     pScreenAdjusted = TRUE,
                     tx2gene = tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.1)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})

#add gene symbols back in
id8 <- match(dex.padj$geneID, gene$ensembl_gene_id)
dex.padj$gene_symbol <- gene$external_gene_name[id8]
dex.padj = dex.padj %>% relocate(gene_symbol, .after = txID)
head(dex.padj)

#export DEXseq_StageR sig. genes/transcripts file
write.csv(dex.padj, file="ccRCC_DEXseq_StageR_padj.csv")
