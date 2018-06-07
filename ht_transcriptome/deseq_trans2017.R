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
library(gridExtra)
library(grid)
library(lattice)

setwd("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/deseq")

# Read in transcript counts and metadata file
#Sample_ids not in same order in the two files. 
trans_cts_lb <- as.matrix(read.csv("~/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/deseq/transcript_count_matrix.csv",
                      row.names="transcript_id"))
coldata_lb <- read.csv("~/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/deseq/Full_pheno_data_2017.csv",
                    header=TRUE, sep=",", row.names="sample_id")

# It is critical that columns of countData and rows of ColData are in the same order AND match! 
all(rownames(coldata_lb) %in% colnames(trans_cts_lb))  #Should return TRUE
countData_lb <- trans_cts_lb[, rownames(coldata_lb)]
all(rownames(coldata_lb) == colnames(countData_lb))    # should return TRUE
#CHOOSE DESIGN 
# Create one Full DESeqDataSet from count matrix and labels
trans_dds_lb <- DESeqDataSetFromMatrix(countData = trans_cts_lb, 
                                       colData = coldata_lb, design = ~ condition)
head(trans_dds_lb)

#PRE-FILTER to remove rows with 0 or 1 read 
trans_dds_lb <- trans_dds_lb[ rowSums(counts(trans_dds_lb)) > 1, ]
head(trans_dds_lb)
#By default, R will choose a reference level for factors based on alphabetical order then DESeq 
#will do the comparisons based on the alphabetical order of the levels.
#RELEVEL: Using relevel, just specifying the reference level
trans_dds_lb$condition <- relevel(trans_dds_lb$condition, ref = "control_1")
#how many genes we capture by counting the number of genes that have non–zero counts in all samples
GeneCounts_trans_lb <- counts(trans_dds_lb)
idx.nz_trans_lb <- apply(GeneCounts_trans_lb, 1, function(x) { all(x > 0)}) 
sum(idx.nz_trans_lb) #5299
#NORMALIZATION: trans
trans_lb <- as.character(colData(trans_dds_lb)$condition)
colData(trans_dds_lb)$condition <- factor(trans_lb, levels = c("control_1", "control_2", "RE_1", "RE_2", 
                                                                "RI_1", "RI_2", "S4_1", "S4_2", "RIplusRE", "S4plusRE"))
trans_dds_lb <- estimateSizeFactors(trans_dds_lb)
sizeFactors(trans_dds_lb)
# C_K_0       C_M_0        C_R1        C_R2        C_R3       C_V_0      RE_K_6      RE_M_6       RE_R1       RE_R2       RE_R3 
# 1.0239149   0.9803377   0.7846750   0.8757910   0.3851861   1.3010847   1.0786567   0.9759637   1.2099079   0.8404709   1.5412151 
# RE_V_6     RI_K_24      RI_K_6     RI_M_24      RI_M_6     RI_V_24      RI_V_6 RIplusRE_R1 RIplusRE_R2 RIplusRE_R3     S4_K_24 
# 2.4005367   1.5422183   1.3627905   1.1836622   1.2417584   1.1618230   1.0777226   0.9634166   1.0389343   0.6546941   1.3536230 
# S4_K_6     S4_M_24      S4_M_6     S4_V_24      S4_V_6 S4plusRE_R1 S4plusRE_R2 S4plusRE_R3 
# 0.7976592   1.0172173   0.8475529   0.6666030   1.0320086   1.1872739   1.4098917   0.8592174 

#Transformation of data
rld_trans_lb <- rlogTransformation(trans_dds_lb, blind=TRUE) 
#rld taking too long: https://support.bioconductor.org/p/77122/
#Suggests using vst instead of rlog transformation for more samples
vst_trans_lb <- varianceStabilizingTransformation(trans_dds_lb, blind=TRUE)
distsRL_trans_lb <- dist(t(assay(vst_trans_lb)))
mat_trans_lb <- as.matrix(distsRL_trans_lb)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat_trans_lb, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "HeatmapPlots_translb.pdf")
dev.off()
#PCA :to visualize sample–to–sample distances
z <- DESeq2::plotPCA(vst_trans_lb, intgroup=c("condition"))
z + geom_label(aes(label = name))

########################################################################################################
#Exp1: The experiment that compares the exposures without pretreatment
########################################################################################################
##prep data for comparison of exp 1: 1,2,6,7,8,12,14,16,18,23,25,27. 6 hour exposures only
#Including all 6h exposures for the analysis then use contrast function in results to get DEGs for each objective. 

trans_cts_1 <- trans_cts_lb[,c(1,2,6,7,8,12,14,16,18,23,25,27)]
coldata_1 <- coldata_lb[c(1,2,6,7,8,12,14,16,18,23,25,27),]
all(rownames(coldata_1) %in% colnames(trans_cts_1))  #Should return TRUE
countData_1 <- trans_cts_1[, rownames(coldata_1)]
all(rownames(coldata_1) == colnames(countData_1))
trans_dds_1 <- DESeqDataSetFromMatrix(countData = trans_cts_1, 
                                       colData = coldata_1, design = ~ condition)
head(trans_dds_1)
trans_dds_1 <- trans_dds_1[ rowSums(counts(trans_dds_1)) > 1, ]
head(trans_dds_1)
trans_dds_1$condition <- relevel(trans_dds_1$condition, ref = "control_1")
GeneCounts_trans_1 <- counts(trans_dds_1)
idx.nz_trans_1 <- apply(GeneCounts_trans_1, 1, function(x) { all(x > 0)}) 
sum(idx.nz_trans_1) #27238
trans_1 <- as.character(colData(trans_dds_1)$condition)
colData(trans_dds_1)$condition <- factor(trans_1, levels = c("control_1", "RE_1", 
                                                               "RI_1", "S4_1"))
trans_dds_1 <- estimateSizeFactors(trans_dds_1)
sizeFactors(trans_dds_1)
# C_K_0     C_M_0     C_V_0    RE_K_6    RE_M_6    RE_V_6    RI_K_6    RI_M_6    RI_V_6    S4_K_6    S4_M_6    S4_V_6 
# 0.9247246 0.9094996 1.1873699 0.9878545 0.9055591 2.2299558 1.2526831 1.1078749 1.0350831 0.7089656 0.7892315 0.9896156 
rld_trans_1 <- rlogTransformation(trans_dds_1, blind=TRUE) 
distsRL_trans_rld_1 <- dist(t(assay(rld_trans_1)))
mat_trans_rld_1 <- as.matrix(distsRL_trans_rld_1)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat_trans_rld_1, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "HeatmapPlots_transrld1.pdf")
dev.off()
DESeq2::plotPCA(rld_trans_1, intgroup=c("condition"))
##MDS plot##
mds_1 <- data.frame(cmdscale(mat_trans_rld_1))
mds_1 <- cbind(mds_1, as.data.frame(colData(rld_trans_1)))
mds_1 <- data.frame(cmdscale(mat_trans_rld_1),eig=TRUE,k=2,x.ret=TRUE)
mds_1 <- cbind(mds_1, as.data.frame(colData(rld_trans_1)))
ggplot(mds_1, aes(X1,X2,color=condition)) + geom_point(size=3)
mds_1
# X1         X2  eig k x.ret condition time sizeFactor
# C_K_0   129.00414 -208.36218 TRUE 2  TRUE control_1    6  0.9247246
# C_M_0   127.54957  226.84164 TRUE 2  TRUE control_1    6  0.9094996
# C_V_0  -308.71186   58.49441 TRUE 2  TRUE control_1    6  1.1873699
# RE_K_6  141.44955 -212.31056 TRUE 2  TRUE      RE_1    6  0.9878545
# RE_M_6  139.62610  156.82668 TRUE 2  TRUE      RE_1    6  0.9055591
# RE_V_6 -504.03246  -70.31705 TRUE 2  TRUE      RE_1    6  2.2299558
# RI_K_6  115.43542 -311.63932 TRUE 2  TRUE      RI_1    6  1.2526831
# RI_M_6  124.03323  247.47239 TRUE 2  TRUE      RI_1    6  1.1078749
# RI_V_6 -152.16795   28.31184 TRUE 2  TRUE      RI_1    6  1.0350831
# S4_K_6  150.59676 -155.40619 TRUE 2  TRUE      S4_1    6  0.7089656
# S4_M_6  129.83220  197.45368 TRUE 2  TRUE      S4_1    6  0.7892315
# S4_V_6  -92.61469   42.63467 TRUE 2  TRUE      S4_1    6  0.9896156
#vst transformation#
vst_trans_1 <- varianceStabilizingTransformation(trans_dds_1, blind=TRUE)
distsRL_trans_1 <- dist(t(assay(vst_trans_1)))
mat_trans_1 <- as.matrix(distsRL_trans_1)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat_trans_1, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "HeatmapPlots_trans_vst_1.pdf")
dev.off()
z1 <- DESeq2::plotPCA(vst_trans_1, intgroup=c("condition"))
z1 + geom_label(aes(label = name))

##All bio reps pulled out separately to plot PCA##
mtr <- trans_cts_lb[,c(2,8,15,16,24,25)]
mcol <- coldata_lb[c(2,8,15,16,24,25),]
mdds <- DESeqDataSetFromMatrix(countData = mtr, colData = mcol, design = ~ condition)
mvst <- varianceStabilizingTransformation(mdds, blind=TRUE)
pcam <- DESeq2::plotPCA(mvst, intgroup=c("condition")) + geom_label(aes(label = name))
pcam <- pcam + theme(text = element_text(size=8))
pcaml <- pcam + theme(legend.position="none")
ktr <- trans_cts_lb[,c(1,7,13,14,22,23)]
kcol <- coldata_lb[c(1,7,13,14,22,23),]
kdds <- DESeqDataSetFromMatrix(countData = ktr, colData = kcol, design = ~ condition)
kvst <- varianceStabilizingTransformation(kdds, blind=TRUE)
pcak <- DESeq2::plotPCA(kvst, intgroup=c("condition")) + geom_label(aes(label = name))
pcak <- pcak + theme(text = element_text(size=8))
pcakl <- pcak + theme(legend.position="none")
vtr <- trans_cts_lb[,c(6,12,17,18,26,27)]
vcol <- coldata_lb[c(6,12,17,18,26,27),]
vdds <- DESeqDataSetFromMatrix(countData = vtr, colData = vcol, design = ~ condition)
vvst <- varianceStabilizingTransformation(vdds, blind=TRUE)
pcav <- DESeq2::plotPCA(vvst, intgroup=c("condition")) + geom_label(aes(label = name))
pcav <- pcav + theme(text = element_text(size=8))
pcavl <- pcav + theme(legend.position="none")
pcavl
pcaml
pcakl
grid.arrange(pcaml, pcakl, pcavl, nrow=2)
##
topVarGenes_trans_1 <- head(order(rowVars(assay(vst_trans_1)), decreasing = TRUE), 50) #picking the 50genes with highest variance across samples
mat_trans_1 <- assay(vst_trans_1)[ topVarGenes_trans_1, ]
mat_trans_1 <- mat_trans_1 - rowMeans(mat_trans_1)
anno_trans_1 <- as.data.frame(colData(vst_trans_1)[, c("condition","time")]) 
pheatmap(mat_trans_1, annotation_col = anno_trans_1)

### DE analysis ####
#Just run the whole deseq function together
trans_dds_1 <- DESeq(trans_dds_1)
#Extracting results
DESeq2Res_trans_1_REvsCon <- results(trans_dds_1, pAdjustMethod = "BH", contrast = c("condition", "RE_1", "control_1")) #BBH=Benjamini Hochberg adjustment
DESeq2Res_trans_1_REvsCon
DESeq2Res_trans_1_REvsCon <- DESeq2Res_trans_1_REvsCon[ !is.na(DESeq2Res_trans_1_REvsCon$padj), ]
sig_trans_1_REvsCon <- DESeq2Res_trans_1_REvsCon[ which(DESeq2Res_trans_1_REvsCon$padj < 0.05 ), ]
nonsig_trans_1_REvsCon <- DESeq2Res_trans_1_REvsCon[ which(DESeq2Res_trans_1_REvsCon$padj > 0.05 ), ]
head( sig_trans_1_REvsCon[ order( sig_trans_1_REvsCon$log2FoldChange ), ] ) #head for strongest downregulation
tail( sig_trans_1_REvsCon[ order( sig_trans_1_REvsCon$log2FoldChange ), ] ) #tail for strongest downregulation
summary(sig_trans_1_REvsCon) #1534
summary(nonsig_trans_1_REvsCon)#58106
write.csv( as.data.frame(sig_trans_1_REvsCon), file="Gene_dfSig_trans_REvsCon_0.05.csv") #1534 genes
write.csv( as.data.frame(nonsig_trans_1_REvsCon), file="Gene_dfNonSig_trans_REvsCon_0.05.csv")

DESeq2Res_trans_1_RIvsCon <- results(trans_dds_1, pAdjustMethod = "BH", contrast = c("condition", "RI_1", "control_1")) #BBH=Benjamini Hochberg adjustment
DESeq2Res_trans_1_RIvsCon
DESeq2Res_trans_1_RIvsCon <- DESeq2Res_trans_1_RIvsCon[ !is.na(DESeq2Res_trans_1_RIvsCon$padj), ]
sig_trans_1_RIvsCon <- DESeq2Res_trans_1_RIvsCon[ which(DESeq2Res_trans_1_RIvsCon$padj < 0.05 ), ]
head( sig_trans_1_RIvsCon[ order( sig_trans_1_RIvsCon$log2FoldChange ), ] ) #head for strongest downregulation
tail( sig_trans_1_RIvsCon[ order( sig_trans_1_RIvsCon$log2FoldChange ), ] ) #tail for strongest downregulation
summary(sig_trans_1_RIvsCon) 
write.csv( as.data.frame(sig_trans_1_RIvsCon), file="Gene_dfSig_trans_RIvsCon_0.05.csv") #1550 genes

DESeq2Res_trans_1_S4vsCon <- results(trans_dds_1, pAdjustMethod = "BH", contrast = c("condition", "S4_1", "control_1")) #BBH=Benjamini Hochberg adjustment
DESeq2Res_trans_1_S4vsCon
DESeq2Res_trans_1_S4vsCon <- DESeq2Res_trans_1_S4vsCon[ !is.na(DESeq2Res_trans_1_S4vsCon$padj), ]
sig_trans_1_S4vsCon <- DESeq2Res_trans_1_S4vsCon[ which(DESeq2Res_trans_1_S4vsCon$padj < 0.05 ), ]
head( sig_trans_1_S4vsCon[ order( sig_trans_1_S4vsCon$log2FoldChange ), ] ) #head for strongest downregulation
tail( sig_trans_1_S4vsCon[ order( sig_trans_1_S4vsCon$log2FoldChange ), ] ) #tail for strongest downregulation
summary(sig_trans_1_S4vsCon) 
write.csv( as.data.frame(sig_trans_1_S4vsCon), file="Gene_dfSig_trans_S4vsCon_0.05.csv") #2269

DESeq2Res_trans_1_S4vsRI <- results(trans_dds_1, pAdjustMethod = "BH", contrast = c("condition", "S4_1", "RI_1"))
sig_trans_1_S4vsRI <- DESeq2Res_trans_1_S4vsRI[ which(DESeq2Res_trans_1_S4vsRI$padj < 0.05 ), ]
head( sig_trans_1_S4vsRI[ order( sig_trans_1_S4vsRI$log2FoldChange ), ] ) #head for strongest downregulation
tail( sig_trans_1_S4vsRI[ order( sig_trans_1_S4vsRI$log2FoldChange ), ] ) #tail for strongest downregulation
summary(sig_trans_1_S4vsRI) 
write.csv( as.data.frame(sig_trans_1_S4vsRI), file="Gene_dfSig_trans_S4vsRI_0.05.csv") #1842

DESeq2Res_trans_1_S4vsRE <- results(trans_dds_1, pAdjustMethod = "BH", contrast = c("condition", "S4_1", "RE_1"))
sig_trans_1_S4vsRE <- DESeq2Res_trans_1_S4vsRE[ which(DESeq2Res_trans_1_S4vsRE$padj < 0.05 ), ]
head( sig_trans_1_S4vsRE[ order( sig_trans_1_S4vsRE$log2FoldChange ), ] ) #head for strongest downregulation
tail( sig_trans_1_S4vsRE[ order( sig_trans_1_S4vsRE$log2FoldChange ), ] ) #tail for strongest downregulation
summary(sig_trans_1_S4vsRE) 
write.csv( as.data.frame(sig_trans_1_S4vsRE), file="Gene_dfSig_trans_S4vsRE_0.05.csv") #1757

DESeq2Res_trans_1_RIvsRE <- results(trans_dds_1, pAdjustMethod = "BH", contrast = c("condition", "RI_1", "RE_1"))
sig_trans_1_RIvsRE <- DESeq2Res_trans_1_RIvsRE[ which(DESeq2Res_trans_1_RIvsRE$padj < 0.05 ), ]
head( sig_trans_1_RIvsRE[ order( sig_trans_1_RIvsRE$log2FoldChange ), ] ) #head for strongest downregulation
tail( sig_trans_1_RIvsRE[ order( sig_trans_1_RIvsRE$log2FoldChange ), ] ) #tail for strongest downregulation
summary(sig_trans_1_RIvsRE) 
write.csv( as.data.frame(sig_trans_1_RIvsRE), file="Gene_dfSig_trans_RIvsRE_0.05.csv") #1229

############################################################################################################################################
##prep data for comparison of exp 1: 1,2,6,13:18,22:27. RI,S4 6and 24h exposures only
trans_cts_2 <- trans_cts_lb[,c(1,2,6,13:18,22:27)]
coldata_2 <- coldata_lb[c(1,2,6,13:18,22:27),]
coldata_2$time <- as.factor(coldata_2$time)
all(rownames(coldata_2) %in% colnames(trans_cts_2))  #Should return TRUE
countData_2 <- trans_cts_2[, rownames(coldata_2)]
all(rownames(coldata_2) == colnames(countData_2))
trans_dds_2 <- DESeqDataSetFromMatrix(countData = trans_cts_2, 
                                      colData = coldata_2, design =  ~ condition)

head(trans_dds_2)
trans_dds_2 <- trans_dds_1[ rowSums(counts(trans_dds_2)) > 1, ]
head(trans_dds_2)
trans_dds_2$condition <- relevel(trans_dds_2$condition, ref = "control_1")
GeneCounts_trans_2 <- counts(trans_dds_2)
idx.nz_trans_2 <- apply(GeneCounts_trans_2, 1, function(x) { all(x > 0)}) 
sum(idx.nz_trans_2) #24753
trans_2 <- as.character(colData(trans_dds_2)$condition)
colData(trans_dds_2)$condition <- factor(trans_2, levels = c("control_1", 
                                                             "RI_1","RI_2", "S4_1","S4_2"))
trans_dds_2 <- estimateSizeFactors(trans_dds_2)
sizeFactors(trans_dds_2)
# C_K_0     C_M_0     C_V_0   RI_K_24    RI_K_6   RI_M_24    RI_M_6   RI_V_24    RI_V_6   S4_K_24    S4_K_6   S4_M_24    S4_M_6   S4_V_24 
# 0.9916253 0.9747851 1.2856738 1.3816964 1.3454485 1.1747307 1.1818767 1.1496264 1.1176171 1.2002987 0.7585019 0.9441510 0.8427420 0.5889437 
# S4_V_6 
# 1.0731633
#rld transformation
rld_trans_2 <- rlogTransformation(trans_dds_2, blind=TRUE) 
distsRL_trans_rld_2 <- dist(t(assay(rld_trans_2)))
mat_trans_rld_2 <- as.matrix(distsRL_trans_rld_2)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat_trans_rld_2, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "HeatmapPlots_transrld2.pdf")
dev.off()
z2 <- DESeq2::plotPCA(rld_trans_2, intgroup=c("condition"))
z2 + geom_label(aes(label = name))
#vst transformation
vst_trans_2 <- varianceStabilizingTransformation(trans_dds_2, blind=TRUE)
distsRL_trans_2 <- dist(t(assay(vst_trans_2)))
mat_trans_2 <- as.matrix(distsRL_trans_2)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat_trans_2, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "HeatmapPlots_trans_vst_2.pdf")
dev.off()
z22 <- DESeq2::plotPCA(vst_trans_2, intgroup=c("condition"))
z22 + geom_label(aes(label = name))
#MDS plot##
mds_2 <- data.frame(cmdscale(mat_trans_rld_2))
ggplot(mds_2, aes(X1,X2)) + geom_point(size=3)
mds_2 <- cbind(mds_2, as.data.frame(colData(rld_trans_2)))
mds_2 <- data.frame(cmdscale(mat_trans_rld_2),eig=TRUE,k=2,x.ret=TRUE)
mds_2 <- cbind(mds_2, as.data.frame(colData(rld_trans_2)))
ggplot(mds_2, aes(X1,X2,color=condition)) + geom_point(size=3)
mds_2
# X1          X2  eig k x.ret condition time sizeFactor
# C_K_0    122.94289  179.170054 TRUE 2  TRUE control_1    6  0.9916253
# C_M_0     54.89742 -200.392543 TRUE 2  TRUE control_1    6  0.9747851
# C_V_0   -446.91990  137.785451 TRUE 2  TRUE control_1    6  1.2856738
# RI_K_24  177.96020  252.683903 TRUE 2  TRUE      RI_2   24  1.3816964
# RI_K_6   127.29169  202.322733 TRUE 2  TRUE      RI_1    6  1.3454485
# RI_M_24   77.80853 -209.244319 TRUE 2  TRUE      RI_2   24  1.1747307
# RI_M_6    60.45312 -212.281346 TRUE 2  TRUE      RI_1    6  1.1818767
# RI_V_24 -210.23913   22.561413 TRUE 2  TRUE      RI_2   24  1.1496264
# RI_V_6  -210.54879    8.709345 TRUE 2  TRUE      RI_1    6  1.1176171
# S4_K_24  174.26016  167.698957 TRUE 2  TRUE      S4_2   24  1.2002987
# S4_K_6   126.55814   84.039540 TRUE 2  TRUE      S4_1    6  0.7585019
# S4_M_24   89.89228 -182.172439 TRUE 2  TRUE      S4_2   24  0.9441510
# S4_M_6    64.11576 -203.829866 TRUE 2  TRUE      S4_1    6  0.8427420
# S4_V_24  -56.78560  -21.387508 TRUE 2  TRUE      S4_2   24  0.5889437
# S4_V_6  -151.68678  -25.663374 TRUE 2  TRUE      S4_1    6  1.0731633
#clustering with vst transformation
topVarGenes_trans_2 <- head(order(rowVars(assay(vst_trans_2)), decreasing = TRUE), 50) #picking the 50genes with highest variance across samples
mat_trans_2 <- assay(vst_trans_2)[ topVarGenes_trans_2, ]
mat_trans_2 <- mat_trans_2 - rowMeans(mat_trans_2)
anno_trans_2 <- as.data.frame(colData(vst_trans_2)[, c("condition","time")]) 
pheatmap(mat_trans_2, annotation_col = anno_trans_2)
#clustering with rld transformation
topVarGenes_rld_trans_2 <- head(order(rowVars(assay(rld_trans_2)), decreasing = TRUE), 50) #picking the 50genes with highest variance across samples
mat_rld_trans_2 <- assay(rld_trans_2)[ topVarGenes_rld_trans_2, ]
mat_rld_trans_2 <- mat_rld_trans_2 - rowMeans(mat_rld_trans_2)
anno_rld_trans_2 <- as.data.frame(colData(rld_trans_2)[, c("condition","time")]) 
pheatmap(mat_rld_trans_2, annotation_col = anno_rld_trans_2)
##
### DE analysis ####
#Just run the whole deseq function together
trans_dds_2 <- DESeq(trans_dds_2)
#Extracting results

DESeq2Res_trans_2_RI_1vsCon <- results(trans_dds_2, pAdjustMethod = "BH", contrast = c("condition", "RI_1", "control_1")) #BBH=Benjamini Hochberg adjustment
DESeq2Res_trans_2_RI_1vsCon
DESeq2Res_trans_2_RI_1vsCon <- DESeq2Res_trans_2_RI_1vsCon[ !is.na(DESeq2Res_trans_2_RI_1vsCon$padj), ]
sig_trans_2_RI_1vsCon <- DESeq2Res_trans_2_RI_1vsCon[ which(DESeq2Res_trans_2_RI_1vsCon$padj < 0.05 ), ]
summary(sig_trans_2_RI_1vsCon) #1806
write.csv( as.data.frame(sig_trans_2_RI_1vsCon), file="Gene_dfSig_trans2_RI_1vsCon_0.05.csv") #1806 genes

DESeq2Res_trans_2_RI_2vsCon <- results(trans_dds_2, pAdjustMethod = "BH", contrast = c("condition", "RI_2", "control_1")) #BBH=Benjamini Hochberg adjustment
DESeq2Res_trans_2_RI_2vsCon
DESeq2Res_trans_2_RI_2vsCon <- DESeq2Res_trans_2_RI_2vsCon[ !is.na(DESeq2Res_trans_2_RI_2vsCon$padj), ]
sig_trans_2_RI_2vsCon <- DESeq2Res_trans_2_RI_2vsCon[ which(DESeq2Res_trans_2_RI_2vsCon$padj < 0.05 ), ]
summary(sig_trans_2_RI_2vsCon) #2139
write.csv( as.data.frame(sig_trans_2_RI_2vsCon), file="Gene_dfSig_trans2_RI_2vsCon_0.05.csv")#2139

DESeq2Res_trans_2_S4_1vsCon <- results(trans_dds_2, pAdjustMethod = "BH", contrast = c("condition", "S4_1", "control_1")) #BBH=Benjamini Hochberg adjustment
DESeq2Res_trans_2_S4_1vsCon
DESeq2Res_trans_2_S4_1vsCon <- DESeq2Res_trans_2_S4_1vsCon[ !is.na(DESeq2Res_trans_2_S4_1vsCon$padj), ]
sig_trans_2_S4_1vsCon <- DESeq2Res_trans_2_S4_1vsCon[ which(DESeq2Res_trans_2_S4_1vsCon$padj < 0.05 ), ]
summary(sig_trans_2_S4_1vsCon) #2413
write.csv( as.data.frame(sig_trans_2_S4_1vsCon), file="Gene_dfSig_trans2_S4_1vsCon_0.05.csv") #2413

DESeq2Res_trans_2_S4_2vsCon <- results(trans_dds_2, pAdjustMethod = "BH", contrast = c("condition", "S4_2", "control_1")) #BBH=Benjamini Hochberg adjustment
DESeq2Res_trans_2_S4_2vsCon
DESeq2Res_trans_2_S4_2vsCon <- DESeq2Res_trans_2_S4_2vsCon[ !is.na(DESeq2Res_trans_2_S4_2vsCon$padj), ]
sig_trans_2_S4_2vsCon <- DESeq2Res_trans_2_S4_2vsCon[ which(DESeq2Res_trans_2_S4_2vsCon$padj < 0.05 ), ]
summary(sig_trans_2_S4_2vsCon) #3459
write.csv( as.data.frame(sig_trans_2_S4_2vsCon), file="Gene_dfSig_trans2_S4_2vsCon_0.05.csv")#3459

########################################################################################################
#Exp2: The experiment that compares the exposures WITH pretreatment
########################################################################################################
##prep data for comparison of exp 1: 3,4,5,9,10,11,19,20,21,28,29,30. 
trans_cts_3 <- trans_cts_lb[,c(3,4,5,9,10,11,19,20,21,28,29,30)]
coldata_3 <- coldata_lb[c(3,4,5,9,10,11,19,20,21,28,29,30),]
head(trans_cts_3)
coldata_3
all(rownames(coldata_3) %in% colnames(trans_cts_3))  #Should return TRUE
countData_3 <- trans_cts_3[, rownames(coldata_3)]
all(rownames(coldata_3) == colnames(countData_3))
trans_dds_3 <- DESeqDataSetFromMatrix(countData = trans_cts_3, 
                                      colData = coldata_3, design = ~ condition)
head(trans_dds_3)
trans_dds_3 <- trans_dds_3[ rowSums(counts(trans_dds_3)) > 1, ]
head(trans_dds_3)
trans_dds_3$condition <- relevel(trans_dds_3$condition, ref = "control_2")
GeneCounts_trans_3 <- counts(trans_dds_3)
idx.nz_trans_3 <- apply(GeneCounts_trans_3, 1, function(x) { all(x > 0)}) 
sum(idx.nz_trans_3) #5982
trans_3 <- as.character(colData(trans_dds_3)$condition)
colData(trans_dds_3)$condition <- factor(trans_3, levels = c("control_2", "RE_2", 
                                                             "RIplusRE", "S4plusRE"))
trans_dds_3 <- estimateSizeFactors(trans_dds_3)
sizeFactors(trans_dds_3)
# C_R1        C_R2        C_R3       RE_R1       RE_R2       RE_R3 RIplusRE_R1 RIplusRE_R2 RIplusRE_R3 S4plusRE_R1 S4plusRE_R2 
# 0.8730196   1.0148181   0.4618319   1.3686746   0.9307247   1.7670194   1.0864780   1.1336485   0.7644821   1.3493908   1.5643641 
# S4plusRE_R3 
# 0.9798007 
rld_trans_3 <- rlogTransformation(trans_dds_3, blind=TRUE) 
distsRL_trans_rld_3 <- dist(t(assay(rld_trans_3)))
mat_trans_rld_3 <- as.matrix(distsRL_trans_rld_3)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat_trans_rld_3, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "HeatmapPlots_transrld_3.pdf")
dev.off()
z3 <- DESeq2::plotPCA(rld_trans_3, intgroup=c("condition"))
z3 + geom_label(aes(label = name))
vst_trans_3 <- varianceStabilizingTransformation(trans_dds_3, blind=TRUE)
distsRL_trans_3 <- dist(t(assay(vst_trans_3)))
mat_trans_3 <- as.matrix(distsRL_trans_3)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat_trans_3, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(pdf, "HeatmapPlots_trans_vst_3.pdf")
dev.off()
z33 <- DESeq2::plotPCA(vst_trans_3, intgroup=c("condition"))
z33 + geom_label(aes(label = name))

topVarGenes_trans_3 <- head(order(rowVars(assay(vst_trans_3)), decreasing = TRUE), 50) #picking the 50genes with highest variance across samples
mat_trans_3 <- assay(vst_trans_3)[ topVarGenes_trans_3, ]
mat_trans_3 <- mat_trans_3 - rowMeans(mat_trans_3)
anno_trans_3 <- as.data.frame(colData(vst_trans_3)[, c("condition","time")]) 
pheatmap(mat_trans_3, annotation_col = anno_trans_3)

