
#Script to generate QC Report of fastq reads
targets <- read.delim("targets.txt")
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/fastqQuality.R")
input1 <- as.character(paste("./data/", targets$FileName1, sep="")); names(input1) <- paste(targets$SampleName,"_","pair1")
input2 <- input2 <- paste("./data/", targets$FileName2, sep=""); names(input2) <- paste(targets$SampleName,"_","pair2")
apt <-append(input1,input2)
        fqlist <- seeFastq(fastq=apt, batchsize=50000, klength=8)
        pdf(paste("./results/","fastqReport.pdf", sep=""), height=18, width=4*length(input1));
        seeFastqPlot(fqlist,cex=1.5); dev.off()
###########################################

#R script to generate Alignment Summary from read alignments
library(ShortRead); library(Rsamtools)
targets <- read.delim("targets.txt")
bfl <- BamFileList(paste0("./results/BAMFiles/", targets$FileName1, ".bam"), yieldSize=50000, index=character())

#For single-end data or if you are creating 2 individual bam files for a pair of sequences(paired-end)
Nreads1 <- ((countLines(paste("./data/", targets$FileName1,sep=""))/4))
Nreads2 <- ((countLines(paste("./data/", targets$FileName2,sep=""))/4))
Nreads = Nreads1 + Nreads2

Nalign <- countBam(bfl)
#For the case where creating 1 single bam file for a pair of sequences(paired end)
(read_statsDF <- data.frame(FileName=names(Nreads), Nreads=Nreads, Nalign=Nalign$records,
                            Perc_Aligned=Nalign$records/Nreads*100))

write.table(read_statsDF,"results/read_statsDF.xls", row.names=FALSE, quote=FALSE, sep="\t")

#Store Annotations in TranscriptDB
library(ShortRead); library(Rsamtools)
library(GenomicFeatures)
txdb <- makeTranscriptDbFromGFF(file="data/TAIR10_GFF3_genes.gff",
        format="gff3",
        dataSource="TAIR",
        species="Arabidopsis thaliana")
saveDb(txdb, file="./data/TAIR10.sqlite")
txdb <- loadDb("./data/TAIR10.sqlite")
eByg <- exonsBy(txdb, by="gene")
print(eByg)

targets <- read.delim("targets.txt")
samples <- as.character(targets$SampleName)
factors_sampl <- as.character(targets$Factor_long)

library(GenomicRanges)
#bfl <- BamFileList(samplespath, yieldSize=50000, index=character())
bfl <- BamFileList(paste0("./results/BAMFiles/", targets$FileName1, ".bam"), yieldSize=50000, index=character())
countDF2 <- summarizeOverlaps(eByg, bfl, mode="Union", ignore.strand=TRUE,singleEnd=FALSE, fragments=FALSE)
countDF2 <- assays(countDF2)$counts
colnames(countDF2) <- factors_sampl
write.table(countDF2, "./results/countDF2.xls", quote=FALSE, sep="\t", col.names = NA)

#Simple RPKM normalization using R
returnRPKM <- function(counts, gffsub) {
        geneLengthsInKB <- sum(width(reduce(gffsub)))/1000 # Length of exon union per gene in kbp
        millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
        rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
        rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
        return(rpkm)
}
countDFrpkm <- apply(countDF2, 2, function(x) returnRPKM(counts=x, gffsub=eByg))
write.table(countDFrpkm, "./results/countDFrpkm.xls", quote=FALSE, sep="\t")

#Generate Sample tree
library(ape)
targets <- read.delim("targets.txt")
#countDFrpkm <- read.delim("./results/countDFrpkm")[,-c(1:2)]
countDFrpkm <- read.delim("./results/countDFrpkm.xls")[,-c(1)]
factors_target <- as.character(targets$Factor_long)
colnames(countDFrpkm) <- factors_target
#countDF <- countDF[, names(samples)]
d <- cor(countDFrpkm, method ="spearman")
hc <- hclust(dist(1-d))
pdf("./results/Sample_tree.pdf")
plot.phylo(as.phylo(hc), type="p",edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
dev.off()

#Filter based on log2 fold change
#countDFrpkm <- cbind(countDFrpkm, log2ratio=log2(countDFrpkm[,1]/countDFrpkm[,2]))
#countDFrpkm <- countDFrpkm[is.finite(countDFrpkm[,3]), ]
#degs2fold <- countDFrpkm[countDFrpkm[,3] >= 1 | countDFrpkm[,3] <= -1,]
#write.table(degs2fold, "./results/degs2fold.xls", quote=FALSE, sep="\t", col.names = NA)
#degs2fold <- read.table("./results/degs2fold.xls")

#Script for making comparisions for differential expression analysis

#labels <- c("Sample3-Sample5","Sample7-Sample13","Sample9-Sample15","Sample11-Sample17", "Sample4-Sample6","Sample8-Sample14","Sample10-Sample16","Sample12-Sample18","Sample2-Sample1","Sample4-Sample3","Sample8-Sample7","Sample10-Sample9", "Sample2-Sample1","Sample6-Sample5","Sample14-Sample13","Sample16-Sample15","Sample18-Sample17","Sample1-Sample3","Sample1-Sample7","Sample1-Sample9","Sample1-Sample11","Sample1-Sample5","Sample1-Sample13","Sample1-Sample15","Sample1-Sample17","Sample2-Sample4","Sample2-Sample8","Sample2-Sample10","Sample2-Sample12","Sample2-Sample6","Sample2-Sample14","Sample2-Sample16","Sample2-Sample18") 
labels <- paste("Sample_", 1:8,"s", sep=""); combn(labels, m=2, FUN=paste, collapse="-")
countDFrpkm <- read.delim("./results/countDFrpkm.xls")[,-c(1)]
countDFrpkm_tmp <- read.delim("./results/countDFrpkm.xls")
colnames(countDFrpkm) <- labels
rownames(countDFrpkm) <- countDFrpkm_tmp[,1]
comp <- c("Sample_1s-Sample_3s","Sample_1s-Sample_5s","Sample_1s-Sample_7s","Sample_2s-Sample_4s","Sample_2s-Sample_6s","Sample_2s-Sample_8s","Sample_2s-Sample_1s","Sample_4s-Sample_3s","Sample_6s-Sample_5s","Sample_8s-Sample_7s")
#comp <- combn(labels, m=2, FUN=paste, collapse="-")
compl <- strsplit(comp, "-")
names(compl) <- comp
r <- sapply(compl, function(x) log2(countDFrpkm[,x[1]]/countDFrpkm[,x[2]]))
rownames(r) <- countDFrpkm_tmp[,1]
deglist <- sapply(colnames(r), function(x) (rownames(r[r[,x] >= 1 | r[,x] <=-1, ])))
deglist <- lapply(deglist, function (x) x[!is.na(x)]) 
#write.table(deglist, "./results/deglist", quote=FALSE, sep="\t")
degfreq <- sapply(rownames(deglist),function(x) freq(x))
degfreq_tmp <- sapply(deglist, length)
deg_count <- as.numeric(degfreq_tmp)
comp_new <- c("T4-5-4_21_D_T4-5-4_38_0.5hr","T4-5-4_21_D_T4-5-4_38_1hr","T4-5-4_21_D_T4-5-4_38_2hr","LL729_21_D_LL729_38_0.5hr","LL729_21_D_LL729_38_1hr","LL729_21_D_LL729_38_2hr","LL729_21_D_T4-5-4_21_D","LL729_38_0.5hr_T4-5-4_38_0.5hr","LL729_38_1hr_T4-5-4_38_1hr","LL729_38_2hr_T4-5-4_38_2hr")

df_degs <- data.frame(Comparisons=comp_new, Counts_2fold = degfreq_tmp)
write.table(df_degs, "./results/deg2fold.xls", quote=FALSE, row.names=FALSE, sep="\t")
pdf("./results/deg_barplot.pdf")
par(las=2)
barplot(deg_count,horiz=TRUE,cex.names=0.8,names.arg=c("1s-3s","1s-5s","1s-7s","2s-4s","2s-6s","2s-8s","2s-1s","4s-3s","6s-5s","8s-7s"))
dev.off()
#df_deglist <- 
#write.table(deglist, "./results/deglist.xls", quote=FALSE, row.names=FALSE, sep="\t")


#################################
## DEG Analysis with glm edgeR ##
#################################
library(edgeR)
targets <- read.delim("targets.txt", comment.char = "#")
targets <- targets[order(targets$Factor),]
table(as.character(targets$Factor))
samples <- as.character(targets$Factor); names(samples) <- paste(as.character(targets$SampleName), "", sep="")
countDF <- read.table("./results/countDF2.xls", check.names=FALSE)
countDF <- countDF[, names(samples)]
countDF[is.na(countDF)] <- 0
group <- as.character(samples)
y <- DGEList(counts=countDF, group=group) # Constructs DGEList object
## Normalization
y <- calcNormFactors(y)
## Contrast matrix is optional but makes analysis more transparent
mycomp <- list(c("T4-5-4_21_D", "T4-5-4_38_0.5hr"), c("T4-5-4_21_D","T4-5-4_38_1hr"), c("T4-5-4_21_D","T4-5-4_38_2hr"), c("LL729_21_D","LL729_38_0.5hr"), c("LL729_21_D","LL729_38_1hr"), c("LL729_21_D","LL729_38_2hr"), c("LL729_21_D","T4-5-4_21_D"), c("LL729_38_0.5hr","T4-5-4_38_0.5hr"), c("LL729_38_1hr","T4-5-4_38_1hr"), c("LL729_38_2hr","T4-5-4_38_2hr"))
myindex <- 1:length(y$samples$group); names(myindex) <- y$samples$group
edgeDF <- data.frame(row.names=rownames(y))
bcv <- 0.1
for(i in seq(along=mycomp)) {
	fit <- exactTest(y[,myindex[mycomp[[i]]]], dispersion=bcv^2)
        deg <- as.data.frame(topTags(fit, n=length(rownames(y))))[,c(1,3,4)]
	colnames(deg) <- paste(paste(mycomp[[i]], collapse="_"), colnames(deg), sep="_")
        edgeDF <- cbind(edgeDF, deg[rownames(y),])
}
## Add functional descriptions
# system("wget ftp://ftp.arabidopsis.org/home/tair/Proteins/TAIR10_functional_descriptions -P ./data/")
desc <- read.delim("/rhome/tgirke/Projects/Julia/RFP-Seq/data/TAIR10_functional_descriptions")[, c(1,3)]
desc[,1] <- gsub("\\..*", "", desc[,1]); desc <- desc[!duplicated(desc[,1]),]
descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
edgeDF <- cbind(Desc=descv[gsub("_.*", "", rownames(edgeDF))], edgeDF)
write.table(edgeDF, "./results/edgeRexact.xls", quote=FALSE, sep="\t", col.names = NA)

## Import
#library(edgeR)
#targets <- read.delim("./data/targets.txt", comment.char = "#")
#targets <- targets[order(targets$SampleName),]
#samples <- as.character(targets$SampleName)
#countDF <- read.delim("./results/countDF2", row.names=1)
#countDF <- countDF[, names(samples)]
#countDF[is.na(countDF)] <- 0
#group <- as.character(samples)
#bcv <- 0.1
#edgeDF <- data.frame(row.names=rownames(countDF))
#y <- DGEList(counts=countDF, group = targets$Factor)
#et <- exactTest(y, dispersion=bcv^2)

#########################
## GO Analysis of DEGs ##
#########################
edgeDF <- read.delim("./results/edgeRexact.xls", row.names=1, check.names=FALSE)
rownames(edgeDF) <- gsub("_.*", "", rownames(edgeDF))
pval <- edgeDF[, grep("_FDR", colnames(edgeDF))]
fold <- edgeDF[, grep("_logFC", colnames(edgeDF))]
pf <- pval <= 0.05 & (fold >= 1 | fold <= -1)
colnames(pf) <- gsub("_FDR", "", colnames(pf))
DEGlist <- sapply(colnames(pf), function(x) rownames(pf[pf[,x],]))
#DEGlist_df <- sapply(DEGlist[[1]][1], function(x) edgeDF[x,])

#write.table(DEGlist_df, "./results/test_DEGlist_df", quote=FALSE, sep="\t", col.names=NA)

library(ggplot2)
df <- data.frame(Comparisons=names(DEGlist), Counts_edgeR_pval_0.05_2fold=sapply(DEGlist, length), Group=2)
write.table(df[,-3], file="results/DEGcounts.xls", quote=FALSE, row.names=FALSE, sep="\t")
pdf("results/DEGcounts.pdf")
ggplot(df, aes(Comparisons, Counts, fill = Group)) + geom_bar(position="dodge") + coord_flip() + opts(axis.text.y=theme_text(angle=0, hjust=1)) + opts(legend.position = "none")
dev.off()
DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
DEGdf <- data.frame(geneID=unlist(DEGlist), CLID=rep(names(DEGlist), sapply(DEGlist, length)), ClusterSize=rep(sapply(DEGlist,length), sapply(DEGlist, length)))
DEGdf_clustmat <- data.frame(geneID=unlist(DEGlist))
#RPKM_val=(edgeDF[geneID,]$colnames(edgeDF)))

#######################################
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/GOHyperGAll.txt")
mydir <- getwd()
setwd("/rhome/tgirke/Projects/Walling/RNA-Seq/data/GO/")
loadData(); load(file="MF_node_affy_list"); load(file="BP_node_affy_list"); load(file="CC_node_affy_list")
setwd(mydir)
BatchResult <- GOCluster_Report(CL_DF=DEGdf, method="all", id_type="gene", CLSZ=10, cutoff=0.01, gocats=c("MF", "BP", "CC"), recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575"))
write.table(BatchResult, "./results/GOBatchResultedgeR.xls", quote=FALSE, sep="\t", col.names = NA)
BatchResultslim <- GOCluster_Report(CL_DF=DEGdf, method="slim", id_type="gene", CLSZ=10, cutoff=0.01, gocats=c("MF", "BP", "CC"), recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575"))
write.table(BatchResultslim, "./results/GOslimBatchResultedgeR.xls", quote=FALSE, sep="\t", col.names = NA)

####################################################################
edgeDF <- read.delim("./results/edgeRexact.xls", row.names=1, check.names=FALSE)
rownames(edgeDF) <- gsub("_.*", "", rownames(edgeDF))
pval <- edgeDF[, grep("_FDR", colnames(edgeDF))]
fold <- edgeDF[, grep("_logFC", colnames(edgeDF))] 
pf <- (fold >= 1 | fold <= -1)
colnames(pf) <- gsub("_FDR", "", colnames(pf))
DEGlist_edgeDF_2fold <- sapply(colnames(pf), function(x) rownames(pf[pf[,x],]))
df_new <- data.frame(Comparisons=names(DEGlist_edgeDF_2fold), Counts_edgeR_2fold=sapply(DEGlist_edgeDF_2fold, length))
write.table(df_new, "./results/DEGlist_edgeDF_2fold.xls", quote=FALSE, row.names=FALSE, sep="\t")

#####################################################################
gene_freq <- cbind(Counts_2fold=df_degs[,2], Counts_edgeR_2fold=df_new[,2], Counts_edgeR_pval_0.05_2fold=df[,2])
df_final <- data.frame(Comparisons=names(DEGlist_edgeDF_2fold),gene_freq)
write.table(df_final, "./results/Counts_final", quote=FALSE, sep="\t",row.names = FALSE)


#######################################################################
edgeDF <- read.delim("./results/edgeRexact.xls", row.names=1, check.names=FALSE)
rownames(edgeDF) <- gsub("_.*", "", rownames(edgeDF))
pval <- edgeDF[, grep("_FDR", colnames(edgeDF))]
fold <- edgeDF[, grep("_logFC", colnames(edgeDF))]
pf <- pval <= 0.05 & (fold >= 4 | fold <= -4)
colnames(pf) <- gsub("_FDR", "", colnames(pf))
DEGlist_edgeDF_4fold <- sapply(colnames(pf), function(x) rownames(pf[pf[,x],]))
df_new_4fold <- data.frame(Comparisons=names(DEGlist_edgeDF_4fold), Counts_edgeR_4fold=sapply(DEGlist_edgeDF_4fold,length))
write.table(df_new_4fold, "./results/DEGlist_edgeDF_4fold.xls", quote=FALSE, row.names=FALSE, sep="\t")

########################################################################
edgeDF <- read.delim("./results/edgeRexact.xls", row.names=1, check.names=FALSE)
rownames(edgeDF) <- gsub("_.*", "", rownames(edgeDF))
pval <- edgeDF[, grep("_FDR", colnames(edgeDF))]
fold <- edgeDF[, grep("_logFC", colnames(edgeDF))]
pf <- pval <= 0.01 & (fold >= 1 | fold <= -1)
colnames(pf) <- gsub("_FDR", "", colnames(pf))
DEGlist_edgeDF_pval_0.01_2fold <- sapply(colnames(pf), function(x) rownames(pf[pf[,x],]))
df_new_pval_0.01_2fold <- data.frame(Comparisons=names(DEGlist_edgeDF_pval_0.01_2fold), Counts_edgeR_pval_0.01_2fold=sapply(DEGlist_edgeDF_pval_0.01_2fold,length))
write.table(df_new_pval_0.01_2fold, "./results/DEGlist_edgeDF_pval_0.01_2fold.xls", quote=FALSE, row.names=FALSE, sep="\t")

#########################################################################
gene_freq <- cbind(Counts_2fold=df_degs[,2], Counts_edgeR_2fold=df_new[,2], Counts_edgeR_pval_0.05_2fold=df[,2], Counts_edgeR_pval_0.01_2fold=df_new_pval_0.01_2fold[,2])
df_final <- data.frame(Comparisons=names(DEGlist_edgeDF_2fold),gene_freq)
write.table(df_final, "./results/Counts_final", quote=FALSE, sep="\t",row.names = FALSE)


edgeDF <- read.delim("./results/edgeRexact.xls", row.names=1, check.names=FALSE)[,-c(1)]
gene_cluster_list <- as.character(unique(DEGdf_clustmat[,1]))
clust_mat <- numeric(0)
#clust_mat <- c()
#gene_list_new <- c("")
for (i in 1:length(gene_cluster_list)) {
	clust_mat <- rbind(clust_mat,edgeDF[gene_cluster_list[i],])
	#gene_list_new <- append(gene_list_new,gene_cluster_list[i])
	#print(clust_mat)
	#print(i)
}
clust_mat_final <- cbind(gene_cluster_list,clust_mat)
clust_df <- data.frame(clust_mat_final)
#row.names(clust_df) <- gene_cluster_list
clust_df_final <- cbind(as.character(clust_df[,1]),clust_df[,2],clust_df[,5],clust_df[,8],clust_df[,11],clust_df[,14],clust_df[,17],clust_df[,20],clust_df[,23],clust_df[,26],clust_df[,29])
clust_df_final_df <- data.frame(clust_df_final)
names(clust_df_final_df) <- c("geneID","T4.5.4_21_D_T4.5.4_38_0.5hr","T4.5.4_21_D_T4.5.4_38_1hr","T4.5.4_21_D_T4.5.4_38_2hr","LL729_21_D_LL729_38_0.5hr","LL729_21_D_LL729_38_1hr","LL729_21_D_LL729_38_2hr","LL729_21_D_T4.5.4_21_D","LL729_38_0.5hr_T4.5.4_38_0.5hr","LL729_38_1hr_T4.5.4_38_1hr","LL729_38_2hr_T4.5.4_38_2hr")
write.table(clust_df_final_df, "./results/Clust_final.xls", quote=FALSE, sep="\t",row.names = FALSE)

#row.names(clust_df) <- gene_cluster_list

#mycomp <- c("SDMSO-S01C1", "SDMSO-S01C4", "SDMSO-S3C1", "SDMSO-S5C4", "RDMSO-R01C1", "RDMSO-R01C4", "RDMSO-R3C1", "RDMSO-R5C4")
#contrasts <- makeContrasts(contrasts=mycomp, levels=design)
#edgeDF <- data.frame(row.names=rownames(y))
#for(i in seq(along=mycomp)) {
	#lrt <- glmLRT(fit, contrast=contrasts[,i]) # Takes DGEGLM object and carries out the likelihood ratio test. 
        #deg <- as.data.frame(topTags(lrt, n=length(rownames(y))))[,c(1,5)]
        #colnames(deg) <- paste(paste(mycomp[i], collapse="_"), colnames(deg), sep="_")
        #edgeDF <- cbind(edgeDF, deg[rownames(y),])
#}

#Differential exon usage analysis using DEXSeq
#source("/rhome/nkatiyar/RNAseq_pipeline_for_Neerja/RNA_seq_workshop_folder/Rrnaseq/data/Fct/gffexonDEXSeq.R")
#source("data/Fct/gffexonDEXSeq.R")
#gff <- import.gff("./data/TAIR10_GFF3_genes.gff", asRangedData=FALSE)
#seqlengths(gff) <- end(ranges(gff[which(elementMetadata(gff)[,"type"]=="chromosome"),]))
#gffexonDEXSeq <- exons2DEXSeq(gff=gff)
#ids <- as.character(elementMetadata(gffexonDEXSeq)[, "ids"])
#countDFdex <- data.frame(row.names=ids)
#samples <- as.character(targets$FileName1)
#samplespath <- paste("./results/BAMFiles/", samples, ".bam", sep="")
#names(samplespath) <- samples
#for(i in samplespath) {
#        aligns <- readBamGappedAlignments(i) # Substitute next two lines with this one.
#        counts <- countOverlaps(gffexonDEXSeq, aligns)
#        countDFdex <- cbind(countDFdex, counts)
#}
#colnames(countDFdex) <- samples
#countDFdex[1:4,1:2]
#write.table(countDFdex, "./results/countDFdex", quote=FALSE, sep="\t", col.names = NA)
#countDFdex <- read.table("./results/countDFdex")

#samples <- as.character(targets$FileName1)
#samplespath <- paste("./results/BAMFiles/", samples, ".bam", sep="")
#names(samplespath) <- samples

