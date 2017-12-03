#Install and load required libraries 

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")
#install.packages("dplyr")
library("dplyr")
#install.packages("fdrtool")
library(fdrtool)
#install.packages("geneplotter")
library("geneplotter")
#install.packages("gplots")
library(gplots)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("genefilter")
library(genefilter)
#install.packages("MDplot")
library(MDplot)
#install.packages("pheatmap")
library(pheatmap)
library(ggplot2)

# Construct Full_PHENO_DATA.csv that contains metadata 
cts <- as.matrix(read.csv("~/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/gene_count_matrix_names.csv", row.names="gene_id"))
coldata <- read.csv("~/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/Full_PHENO_DATA_C_VIR.csv", header=TRUE, sep=",", row.names=1)

# It is critical that columns of countData and rows of ColData are in the same or order AND match! 
#%in% is an operator that checks if elements in one vector are present in another.but this will 
#return TRUE or FALSE for all the elements
#But if you use 'all' that will check if all elements in one vector are in another and return just one TRUE or FALSE
all(rownames(coldata) %in% colnames(cts))  #Should return TRUE
countData <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(countData))    # should return TRUE

#CHOOSE DESIGN 
# Create one Full DESeqDataSet from count matrix and labels 
dds <- DESeqDataSetFromMatrix(countData = cts, 
                                            colData = coldata, design = ~ condition)
head(dds) #why are my domensions only 6x6? 6colnmaes for 6sras but why just 6 rows?

#Pre-filter to remove rows with 0 or 1 read
dds <- dds[ rowSums(counts(dds)) > 1, ]
head(dds)


#By default, R will choose a reference level for factors based on alphabetical order then DESeq 
#will do the comparisons based on the alphabetical order of the levels.
#Using relevel, just specifying the reference level:
dds$condition <- relevel(dds$condition, ref = "control")

#how many genes we capture by counting the number of genes that have non–zero counts in all samples.
GeneCounts <- counts(dds)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

#normalization per tutorial (DESeq2-analysis.pdf)
#the ratios of the size factors should roughly match the ratios of the library sizes. 
#Dividing each column of the count table by the corresponding size factor yields normalized count values, 
#which can be scaled to give a counts per million interpretation
#if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
con <- as.character(colData(dds)$condition)
colData(dds)$condition <- factor(con, levels = c("control", "rif"))
dds <- estimateSizeFactors(dds)
sizeFactors(dds) #sizefactors for all are around 1 so they have been sequenced equally deeply

# plot densities of counts for the different samples to assess their distributions
#plotting only for genes with non zero counts hence using [idx.nz ,] as done on line 47.
multiecdf( counts(dds, normalized = T)[idx.nz ,], xlab="mean counts", xlim=c(0, 1000)) 
multidensity( counts(dds, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
# All samples dont exactly match in the counts and density.. problem?

## check pairwise MA plots
MA.idx = t(combn(1:8, 2)) 
  for( i in 1:15){
          MDPlot(counts(dds, normalized = T)[idx.nz ,],
          c(MA.idx[i,1],MA.idx[i,2]),
          main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
          colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3) ) 
    }
dev.copy(pdf, "pairwiseMAs.pdf")
dev.off()
# no shift seen so good

##PCA and sample heatmaps
#Which samples are similar to each other, which are different? 
#Does this fit to the expectation from the experiment’s design?
# produce rlog-transformed data
#The aim of the regularized log–transform is to stabilize the variance of the data
#and to make its distribution roughly symmetric 
#For genes with high counts, the rlog transformation differs not much from an ordinary log2 transformation.
#For genes with lower counts, however, the values are shrunken towards the genes’ averages across all samples.
rld <- rlogTransformation(dds, blind=TRUE)
## create a distance matrix between the samples
#pdf("HeatmapPlots.pdf")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "HeatmapPlots.pdf")
dev.off()
#PCA :to visualize sample–to–sample distances
DESeq2::plotPCA(rld, intgroup=c("condition"))
#pca analysis shows high variability among the samples. Expected. 
#26,19,23 cluster together, 20 in far right corner, 24 and 25 in the middle (in order from left) 
#19,24:day5, 20,25:day12 and 23,26:day16

#Gene clustering#
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50) #picking the 50genes with highest variance across samples
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("condition","stressorLevel")]) 
pheatmap(mat, annotation_col = anno)

### DE analysis ####
# We start with raw counts again not log transformed as done for PCA.
#estimate dispersion
#the first step in the analysis of differential expression, 
#is to obtain an estimate of the dispersion parameter for each gene. 
#The typical shape of the dispersion fit is an exponentially decaying curve.
dds <- estimateDispersions(dds)
plotDispEsts(dds)
#The black points are the dispersion estimates for each gene as obtained 
#by considering the information from each gene separately.
#Unless one has many samples, these values fluctuate strongly around their true values.
#the red trend line is fitted, which shows the dispersions’ dependence on the mean
#final estimates (blue points) that are then used in the hypothesis test
#The blue circles above the main “cloud” of points are genes which have high gene–wise dispersion estimates 
#which are labelled as dispersion outliers. These estimates are therefore not shrunk toward the fitted trend line.
#The warnings just indicate that the dispersion estimation failed for some genes
#0.01 means expression tends to differ 10% between sampels of same treatment grps.  
#In our case the red line stays around 1e+00. LESS DISPERSION?

###Statistical testing of differential expression
dds <-  nbinomWaldTest(dds)
DESeq2Res <- results(dds, pAdjustMethod = "BH") #BBH=Benjamini Hochberg adjustment
DESeq2Res #DataFrame with 40905 rows and 6 columns
#Here lfcSE: the standard error estimate for the log2 fold change estimate.
#Null hypo:that there is zero effect of the treatment on the gene and that the observed difference between 
#treatment and control was merely caused by experimental variability.
#pvalue:p value indicates the probability that a fold change as strong as the observed one, or even stronger, 
#would be seen under the situation described by the null hypothesis.
head(DESeq2Res)
summary(DESeq2Res)
#out of 40905 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 32, 0.078% 
#LFC < 0 (down)   : 59, 0.14% 
#outliers [1]     : 5368, 13% 
#low counts [2]   : 793, 1.9% 
#(mean count < 1)
### number of siginificant DE-genes
table(DESeq2Res$padj < 0.1) #identified 91 differentially expressed genes at padj > 0.1
table(DESeq2Res$padj < 0.05) #identified 58 differentially expressed genes at padj < 0.05

# ## get average expressions
# overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])
# 
# ##independent filtering
# attr(DESeq2Res,"filterThreshold") #returns NULL!
# plot(attr(DESeq2Res,"filterNumRej"),type="b", xlab="quantiles of 'baseMean'",
#      ylab="number of rejections") #gives warnings probably coz value of lin e 129 is NULL!

#Inspection and correction of p–values
hist(DESeq2Res$pvalue, col = "lavender", main = "control vs rif", xlab = "p-values")

### remove filtered out genes by independent filtering, ### they have NA adj. pvals
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]

### remove genes with NA pvals (outliers)
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]
### remove adjsuted pvalues, since we add the fdrtool results later on ### (based on the correct p-values)
DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
### use z-scores as input to FDRtool to re-estimate the p-value
FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)

FDR.DESeq2Res$param[1, "sd"] #sd is 0.9357207
DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
hist(FDR.DESeq2Res$pval, col = "royalblue4",
    main = "control vs rif, correct null model", xlab = "CORRECTED p-values")
table(DESeq2Res[,"padj"] < 0.1) #146 differentially expressed genes at a FDR of 0.1. 
table(DESeq2Res[,"padj"] < 0.05) #102 differentially expressed genes at a FDR of 0.1. 
plotMA(DESeq2Res)

#Subset the results table to the differentially expressed genes under FDR 0.1, 
#order the Log2FC table first by strongest down regulation
sig <- DESeq2Res[ which(DESeq2Res$padj < 0.05 ), ]
head( sig[ order( sig$log2FoldChange ), ] ) #head for strongest downregulation
tail( sig[ order( sig$log2FoldChange ), ] ) #tail for strongest downregulation
summary(sig)
nonsig <- DESeq2Res[ which(DESeq2Res$padj > 0.05 ), ]
summary(nonsig)

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(sig)
plotMA(nonsig)

#Export Results to CSV
write.csv( as.data.frame(DESeq2Res), file="Gene_df2.csv")
write.csv( as.data.frame(sig), file="Gene_dfSig2.csv")
write.csv( as.data.frame(nonsig), file = "Gene_df_non_Sig2.csv")

#Session info for records. 
devtools::session_info()

#################################################################################################################

#CHOOSE DESIGN x: Time/StressorLevel  
# Create one Full DESeqDataSet from count matrix and labels 
dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = coldata, design = ~ stressorLevel)
head(dds) #why are my domensions only 6x6? 6colnmaes for 6sras but why just 6 rows?

#Pre-filter to remove rows with 0 or 1 read
dds <- dds[ rowSums(counts(dds)) > 1, ]
head(dds)

#By default, R will choose a reference level for factors based on alphabetical order then DESeq 
#will do the comparisons based on the alphabetical order of the levels.
#Using relevel, just specifying the reference level:
dds$condition <- relevel(dds$condition, ref = "control")

#how many genes we capture by counting the number of genes that have non–zero counts in all samples.
GeneCounts <- counts(dds)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

#normalization per tutorial (DESeq2-analysis.pdf)
#the ratios of the size factors should roughly match the ratios of the library sizes. 
#Dividing each column of the count table by the corresponding size factor yields normalized count values, 
#which can be scaled to give a counts per million interpretation
#if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
con <- as.character(colData(dds)$stressorLevel)
colData(dds)$stressorLevel <- factor(con, levels = c("5d", "12d","16d"))
dds <- estimateSizeFactors(dds)
sizeFactors(dds) #sizefactors for all are around 1 so they have been sequenced equally deeply

# plot densities of counts for the different samples to assess their distributions
#plotting only for genes with non zero counts hence using [idx.nz ,] as done on line 47.
multiecdf( counts(dds, normalized = T)[idx.nz ,], xlab="mean counts", xlim=c(0, 1000)) 
multidensity( counts(dds, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
# All samples dont exactly match in the counts and density.. problem?

## check pairwise MA plots
MA.idx = t(combn(1:8, 2)) 
for( i in 1:15){
  MDPlot(counts(dds, normalized = T)[idx.nz ,],
         c(MA.idx[i,1],MA.idx[i,2]),
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3) ) 
}
dev.copy(pdf, "pairwiseMAs.pdf")
dev.off()
# no shift seen so good

##PCA and sample heatmaps
#Which samples are similar to each other, which are different? 
#Does this fit to the expectation from the experiment’s design?
# produce rlog-transformed data
#The aim of the regularized log–transform is to stabilize the variance of the data
#and to make its distribution roughly symmetric 
#For genes with high counts, the rlog transformation differs not much from an ordinary log2 transformation.
#For genes with lower counts, however, the values are shrunken towards the genes’ averages across all samples.
rld <- rlogTransformation(dds, blind=TRUE)
## create a distance matrix between the samples
#pdf("HeatmapPlots.pdf")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "Design2HeatmapPlots.pdf")
dev.off()
#PCA :to visualize sample–to–sample distances
DESeq2::plotPCA(rld, intgroup=c("condition"))
#pca analysis shows high variability among the samples. Expected. 
#26,19,23 cluster together, 20 in far right corner, 24 and 25 in the middle (in order from left) 
#19,24:day5, 20,25:day12 and 23,26:day16

#estimate dispersion
#the first step in the analysis of differential expression, 
#is to obtain an estimate of the dispersion parameter for each gene. 
#The typical shape of the dispersion fit is an exponentially decaying curve.
dds <- estimateDispersions(dds)
plotDispEsts(dds)
#The black points are the dispersion estimates for each gene as obtained 
#by considering the information from each gene separately.
#Unless one has many samples, these values fluctuate strongly around their true values.
#the red trend line is fitted, which shows the dispersions’ dependence on the mean
#final estimates (blue points) that are then used in the hypothesis test
#The blue circles above the main “cloud” of points are genes which have high gene–wise dispersion estimates 
#which are labelled as dispersion outliers. These estimates are therefore not shrunk toward the fitted trend line.
#The warnings just indicate that the dispersion estimation failed for some genes
#0.01 means expression tends to differ 10% between sampels of same treatment grps.  
#In our case the red line stays around 1e+00. LESS DISPERSION?

###Statistical testing of differential expression
dds <-  nbinomWaldTest(dds)
DESeq2Res <- results(dds, pAdjustMethod = "BH")
head(DESeq2Res)
### number of siginificant DE-genes
table(DESeq2Res$padj < 0.1) #identified 1418 differentially expressed genes

# ## get average expressions
# overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])
# 
# ##independent filtering
# attr(DESeq2Res,"filterThreshold") #returns NULL!
# plot(attr(DESeq2Res,"filterNumRej"),type="b", xlab="quantiles of 'baseMean'",
#      ylab="number of rejections") #gives warnings probably coz value of lin e 129 is NULL!

#Inspection and correction of p–values
hist(DESeq2Res$pvalue, col = "lavender", main = "Time of exposure", xlab = "p-values")

### remove filtered out genes by independent filtering, ### they have NA adj. pvals
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]

### remove genes with NA pvals (outliers)
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]
### remove adjsuted pvalues, since we add the fdrtool results later on ### (based on the correct p-values)
DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
### use z-scores as input to FDRtool to re-estimate the p-value
FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)

FDR.DESeq2Res$param[1, "sd"] #sd is 0.9357
DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
hist(FDR.DESeq2Res$pval, col = "royalblue4",
     main = "control vs rif, correct null model", xlab = "CORRECTED p-values")
table(DESeq2Res[,"padj"] < 0.1) #1874 differentially expressed genes at a FDR of 0.1. 
plotMA(DESeq2Res)

#Subset the results table to the differentially expressed genes under FDR 0.1, 
#order the Log2FC table first by strongest down regulation
sig <- DESeq2Res[ which(DESeq2Res$padj < 0.05 ), ]
head( sig[ order( sig$log2FoldChange ), ] ) #head for strongest downregulation
tail( sig[ order( sig$log2FoldChange ), ] ) #tail for strongest downregulation
summary(sig)
nonsig <- DESeq2Res[ which(DESeq2Res$padj > 0.05 ), ]
summary(nonsig)

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(sig)
plotMA(nonsig)

###################################################################################
#Design2: remove con 12d and run the same analysis as design1 
###################################################################################
# Full_PHENO_DATA_des2.csv that contains metadata with Con 12d removed.
cts <- as.matrix(read.csv("~/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design2/gene_count_matrix_names_des2.csv", row.names="gene_id"))
coldata <- read.csv("~/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design2/Full_PHENO_DATA_C_VIR_des2.csv", header=TRUE, sep=",", row.names=1)

# It is critical that columns of countData and rows of ColData are in the same or order AND match! 
#%in% is an operator that checks if elements in one vector are present in another.but this will 
#return TRUE or FALSE for all the elements
#But if you use 'all' that will check if all elements in one vector are in another and return just one TRUE or FALSE
all(rownames(coldata) %in% colnames(cts))  #Should return TRUE
countData <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(countData))    # should return TRUE

#CHOOSE DESIGN 
# Create one Full DESeqDataSet from count matrix and labels 
dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = coldata, design = ~ condition)
head(dds) 

#Pre-filter to remove rows with 0 or 1 read
dds <- dds[ rowSums(counts(dds)) > 1, ]
head(dds)


#By default, R will choose a reference level for factors based on alphabetical order then DESeq 
#will do the comparisons based on the alphabetical order of the levels.
#Using relevel, just specifying the reference level:
dds$condition <- relevel(dds$condition, ref = "control")

#how many genes we capture by counting the number of genes that have non–zero counts in all samples.
GeneCounts <- counts(dds)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

#normalization per tutorial (DESeq2-analysis.pdf)
#the ratios of the size factors should roughly match the ratios of the library sizes. 
#Dividing each column of the count table by the corresponding size factor yields normalized count values, 
#which can be scaled to give a counts per million interpretation
#if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
con <- as.character(colData(dds)$condition)
colData(dds)$condition <- factor(con, levels = c("control", "rif"))
dds <- estimateSizeFactors(dds)
sizeFactors(dds) #sizefactors for all are around 1 so they have been sequenced equally deeply?
#con5d    con16d     rif5d    rif12d    rif16d 
#0.9742111 1.1075949 1.7072417 0.7324887 0.7402491

# plot densities of counts for the different samples to assess their distributions
#plotting only for genes with non zero counts hence using [idx.nz ,] as done on line 47.
multiecdf( counts(dds, normalized = T)[idx.nz ,], xlab="mean counts", xlim=c(0, 1000)) 
multidensity( counts(dds, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
# All samples dont exactly match in the counts and density.. problem?

## check pairwise MA plots
MA.idx = t(combn(1:8, 2)) 
for( i in 1:15){
  MDPlot(counts(dds, normalized = T)[idx.nz ,],
         c(MA.idx[i,1],MA.idx[i,2]),
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3) ) 
}
dev.copy(pdf, "pairwiseMAs_des2.pdf")
dev.off()
# no shift seen so good

##PCA and sample heatmaps
#Which samples are similar to each other, which are different? 
#Does this fit to the expectation from the experiment’s design?
# produce rlog-transformed data
#The aim of the regularized log–transform is to stabilize the variance of the data
#and to make its distribution roughly symmetric 
#For genes with high counts, the rlog transformation differs not much from an ordinary log2 transformation.
#For genes with lower counts, however, the values are shrunken towards the genes’ averages across all samples.
rld <- rlogTransformation(dds, blind=TRUE)
## create a distance matrix between the samples
#pdf("HeatmapPlots.pdf")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "HeatmapPlots_des2.pdf")
dev.off()
#PCA :to visualize sample–to–sample distances
DESeq2::plotPCA(rld, intgroup=c("condition"))
#pca analysis shows high variability among the samples. Expected. 
#26,19,23 cluster together, 20 in far right corner, 24 and 25 in the middle (in order from left) 
#19,24:day5, 20,25:day12 and 23,26:day16

#Gene clustering#
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50) #picking the 50genes with highest variance across samples
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("condition","stressorLevel")]) 
pheatmap(mat, annotation_col = anno)

### DE analysis ####
# We start with raw counts again not log transformed as done for PCA.
#estimate dispersion
#the first step in the analysis of differential expression, 
#is to obtain an estimate of the dispersion parameter for each gene. 
#The typical shape of the dispersion fit is an exponentially decaying curve.
dds <- estimateDispersions(dds)
plotDispEsts(dds)
#The black points are the dispersion estimates for each gene as obtained 
#by considering the information from each gene separately.
#Unless one has many samples, these values fluctuate strongly around their true values.
#the red trend line is fitted, which shows the dispersions’ dependence on the mean
#final estimates (blue points) that are then used in the hypothesis test
#The blue circles above the main “cloud” of points are genes which have high gene–wise dispersion estimates 
#which are labelled as dispersion outliers. These estimates are therefore not shrunk toward the fitted trend line.
#The warnings just indicate that the dispersion estimation failed for some genes
#0.01 means expression tends to differ 10% between sampels of same treatment grps.  
#In our case the red line stays around 1e+00. LESS DISPERSION?

###Statistical testing of differential expression
dds <-  nbinomWaldTest(dds)
DESeq2Res <- results(dds, pAdjustMethod = "BH") #BH=Benjamini Hochberg adjustment
DESeq2Res #DataFrame with 40905 rows and 6 columns
#Here lfcSE: the standard error estimate for the log2 fold change estimate.
#Null hypo:that there is zero effect of the treatment on the gene and that the observed difference between 
#treatment and control was merely caused by experimental variability.
#pvalue:p value indicates the probability that a fold change as strong as the observed one, or even stronger, 
#would be seen under the situation described by the null hypothesis.
head(DESeq2Res)
summary(DESeq2Res)
#out of 40905 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 32, 0.078% 
#LFC < 0 (down)   : 59, 0.14% 
#outliers [1]     : 5368, 13% 
#low counts [2]   : 793, 1.9% 
#(mean count < 1)
### number of siginificant DE-genes
table(DESeq2Res$padj < 0.1) #identified 2876 differentially expressed genes at padj > 0.1
table(DESeq2Res$padj < 0.05) #identified 1653 differentially expressed genes at padj < 0.05

# ## get average expressions
# overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])
# 
# ##independent filtering
# attr(DESeq2Res,"filterThreshold") #returns NULL!
# plot(attr(DESeq2Res,"filterNumRej"),type="b", xlab="quantiles of 'baseMean'",
#      ylab="number of rejections") #gives warnings probably coz value of lin e 129 is NULL!

#Inspection and correction of p–values
hist(DESeq2Res$pvalue, col = "lavender", main = "control vs rif", xlab = "p-values")

### remove filtered out genes by independent filtering, ### they have NA adj. pvals
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]

### remove genes with NA pvals (outliers)
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]
### remove adjsuted pvalues, since we add the fdrtool results later on ### (based on the correct p-values)
DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
### use z-scores as input to FDRtool to re-estimate the p-value
FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)

FDR.DESeq2Res$param[1, "sd"] #sd is 1.453482
DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
hist(FDR.DESeq2Res$pval, col = "royalblue4",
     main = "control vs rif, correct null model des2", xlab = "CORRECTED p-values")
table(DESeq2Res[,"padj"] < 0.1) #31 differentially expressed genes at a FDR of 0.1. 
table(DESeq2Res[,"padj"] < 0.05) #30 differentially expressed genes at a FDR of 0.1. 
plotMA(DESeq2Res)

#Subset the results table to the differentially expressed genes under FDR 0.1, 
#order the Log2FC table first by strongest down regulation
sig <- DESeq2Res[ which(DESeq2Res$padj < 0.05 ), ]
head( sig[ order( sig$log2FoldChange ), ] ) #head for strongest downregulation
tail( sig[ order( sig$log2FoldChange ), ] ) #tail for strongest upregulation
summary(sig)
nonsig <- DESeq2Res[ which(DESeq2Res$padj > 0.05 ), ]
summary(nonsig)

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(sig)
plotMA(nonsig)

#Export Results to CSV
write.csv( as.data.frame(DESeq2Res), file="Gene_df_des2.csv")
write.csv( as.data.frame(sig), file="Gene_dfSig2_des2.csv")
write.csv( as.data.frame(nonsig), file = "Gene_df_non_Sig2_des2.csv")

#Session info for records. 
devtools::session_info()

####################################################################################
#CHOOSE DESIGN 3: Time course evaluation 
####################################################################################
# Construct Full_PHENO_DATA.csv that contains metadata 
cts <- as.matrix(read.csv("~/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/gene_count_matrix_names.csv", row.names="gene_id"))
coldata <- read.csv("~/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/Full_PHENO_DATA_C_VIR.csv", header=TRUE, sep=",", row.names=1)
coldata$stressorLevel <- as.factor(coldata$stressorLevel)

# It is critical that columns of countData and rows of ColData are in the same or order AND match! 
#%in% is an operator that checks if elements in one vector are present in another.but this will 
#return TRUE or FALSE for all the elements
#But if you use 'all' that will check if all elements in one vector are in another and return just one TRUE or FALSE
all(rownames(coldata) %in% colnames(cts))  #Should return TRUE
countData <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(countData))    # should return TRUE

# Create one Full DESeqDataSet from count matrix and labels 
ddsTC <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = coldata, design = ~ condition + stressorLevel + condition:stressorLevel)
head(ddsTC) #why are my domensions only 6x6? 6colnmaes for 6sras but why just 6 rows?

ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ condition + stressorLevel)
# Warning message:
#   In checkForExperimentalReplicates(object, modelMatrix) :
#   same number of samples and coefficients to fit,
# estimating dispersion by treating samples as replicates.
# read the ?DESeq section on 'Experiments without replicates'
# Experiments without replicates do not allow for estimation of the dispersion of counts around the expected value for each group, 
#which is critical for differential expression analysis.
#all the samples are considered as replicates of a single group for the estimation of dispersion
# "Some overestimation of the variance may be expected, which will make that approach conservative." Furthermore, "while one may not want to draw strong conclusions from such an analysis, 
# it may still be useful for exploration and hypothesis generation."
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),], 4)
# log2 fold change (MLE): conditionrif.stressorLevel5d 
# LRT p-value: '~ condition + stressorLevel + condition:stressorLevel' vs '~ condition + stressorLevel' 
# DataFrame with 4 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat    pvalue      padj
# <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
#   gene33415 71.772216      -5.422973  3.497668  2.413608 0.2991519 0.9603658
# gene33380  7.687932      -5.997684  3.963931  3.908111 0.1416983 0.9603658
# gene33417  3.768397       8.818171  7.059547  2.329659 0.3119758 0.9603658
# gene33416  8.374238      13.051802  7.397494  2.633301 0.2680316 0.9603658

#plot the counts for the groups over time
#for the gene with the smallest adjusted p value, 
#testing for condition-dependent time profile and accounting for differences at time 0
#Keep in mind that the interaction terms are the difference between the two groups 
#at a given time after accounting for the difference at time 0.
tc <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("stressorLevel","condition"), returnData = TRUE)
ggplot(tc,
       aes(x = as.numeric(stressorLevel), y = count, color = condition, group = condition)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10()
#Getting warning messages for this step. "fewer data values than degrees of freedom" and others. 

#Wald tests for the log2 fold changes at individual time points can be investigated using the test argument to results:
resultsNames(ddsTC)
# [1] "Intercept"                    "condition_rif_vs_control"     "stressorLevel_12_vs_5"        "stressorLevel_16_vs_5"       
# [5] "conditionrif.stressorLevel12" "conditionrif.stressorLevel16"
res30 <- results(ddsTC, name="conditionrif.stressorLevel12", test="Wald")
res30[which.min(resTC$padj),]

#We can furthermore cluster significant genes by their profiles. 
#We extract a matrix of the shrunken log2 fold changes using the coef function:
betas <- coef(ddsTC)
colnames(betas)
topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)


#****#
dds <- dds[ rowSums(counts(dds)) > 1, ]
head(dds)

#By default, R will choose a reference level for factors based on alphabetical order then DESeq 
#will do the comparisons based on the alphabetical order of the levels.
#Using relevel, just specifying the reference level:
dds$condition <- relevel(dds$condition, ref = "control")

#how many genes we capture by counting the number of genes that have non–zero counts in all samples.
GeneCounts <- counts(dds)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

#normalization per tutorial (DESeq2-analysis.pdf)
#the ratios of the size factors should roughly match the ratios of the library sizes. 
#Dividing each column of the count table by the corresponding size factor yields normalized count values, 
#which can be scaled to give a counts per million interpretation
#if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
con <- as.character(colData(dds)$stressorLevel)
colData(dds)$stressorLevel <- factor(con, levels = c("5d", "12d","16d"))
dds <- estimateSizeFactors(dds)
sizeFactors(dds) #sizefactors for all are around 1 so they have been sequenced equally deeply

# plot densities of counts for the different samples to assess their distributions
#plotting only for genes with non zero counts hence using [idx.nz ,] as done on line 47.
multiecdf( counts(dds, normalized = T)[idx.nz ,], xlab="mean counts", xlim=c(0, 1000)) 
multidensity( counts(dds, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
# All samples dont exactly match in the counts and density.. problem?

## check pairwise MA plots
MA.idx = t(combn(1:8, 2)) 
for( i in 1:15){
  MDPlot(counts(dds, normalized = T)[idx.nz ,],
         c(MA.idx[i,1],MA.idx[i,2]),
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3) ) 
}
dev.copy(pdf, "pairwiseMAs.pdf")
dev.off()
# no shift seen so good

##PCA and sample heatmaps
#Which samples are similar to each other, which are different? 
#Does this fit to the expectation from the experiment’s design?
# produce rlog-transformed data
#The aim of the regularized log–transform is to stabilize the variance of the data
#and to make its distribution roughly symmetric 
#For genes with high counts, the rlog transformation differs not much from an ordinary log2 transformation.
#For genes with lower counts, however, the values are shrunken towards the genes’ averages across all samples.
rld <- rlogTransformation(dds, blind=TRUE)
## create a distance matrix between the samples
#pdf("HeatmapPlots.pdf")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "Design2HeatmapPlots.pdf")
dev.off()
#PCA :to visualize sample–to–sample distances
DESeq2::plotPCA(rld, intgroup=c("condition"))
#pca analysis shows high variability among the samples. Expected. 
#26,19,23 cluster together, 20 in far right corner, 24 and 25 in the middle (in order from left) 
#19,24:day5, 20,25:day12 and 23,26:day16

#estimate dispersion
#the first step in the analysis of differential expression, 
#is to obtain an estimate of the dispersion parameter for each gene. 
#The typical shape of the dispersion fit is an exponentially decaying curve.
dds <- estimateDispersions(dds)
plotDispEsts(dds)
#The black points are the dispersion estimates for each gene as obtained 
#by considering the information from each gene separately.
#Unless one has many samples, these values fluctuate strongly around their true values.
#the red trend line is fitted, which shows the dispersions’ dependence on the mean
#final estimates (blue points) that are then used in the hypothesis test
#The blue circles above the main “cloud” of points are genes which have high gene–wise dispersion estimates 
#which are labelled as dispersion outliers. These estimates are therefore not shrunk toward the fitted trend line.
#The warnings just indicate that the dispersion estimation failed for some genes
#0.01 means expression tends to differ 10% between sampels of same treatment grps.  
#In our case the red line stays around 1e+00. LESS DISPERSION?

###Statistical testing of differential expression
dds <-  nbinomWaldTest(dds)
DESeq2Res <- results(dds, pAdjustMethod = "BH")
head(DESeq2Res)
### number of siginificant DE-genes
table(DESeq2Res$padj < 0.1) #identified 1418 differentially expressed genes

# ## get average expressions
# overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])
# 
# ##independent filtering
# attr(DESeq2Res,"filterThreshold") #returns NULL!
# plot(attr(DESeq2Res,"filterNumRej"),type="b", xlab="quantiles of 'baseMean'",
#      ylab="number of rejections") #gives warnings probably coz value of lin e 129 is NULL!

#Inspection and correction of p–values
hist(DESeq2Res$pvalue, col = "lavender", main = "Time of exposure", xlab = "p-values")

### remove filtered out genes by independent filtering, ### they have NA adj. pvals
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]

### remove genes with NA pvals (outliers)
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]
### remove adjsuted pvalues, since we add the fdrtool results later on ### (based on the correct p-values)
DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
### use z-scores as input to FDRtool to re-estimate the p-value
FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)

FDR.DESeq2Res$param[1, "sd"] #sd is 0.9357
DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
hist(FDR.DESeq2Res$pval, col = "royalblue4",
     main = "control vs rif, correct null model", xlab = "CORRECTED p-values")
table(DESeq2Res[,"padj"] < 0.1) #1874 differentially expressed genes at a FDR of 0.1. 
plotMA(DESeq2Res)

#Subset the results table to the differentially expressed genes under FDR 0.1, 
#order the Log2FC table first by strongest down regulation
sig <- DESeq2Res[ which(DESeq2Res$padj < 0.05 ), ]
head( sig[ order( sig$log2FoldChange ), ] ) #head for strongest downregulation
tail( sig[ order( sig$log2FoldChange ), ] ) #tail for strongest downregulation
summary(sig)
nonsig <- DESeq2Res[ which(DESeq2Res$padj > 0.05 ), ]
summary(nonsig)

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(sig)
plotMA(nonsig)

####################################################################################

