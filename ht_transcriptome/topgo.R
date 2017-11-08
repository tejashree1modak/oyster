#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")
library("topGO")
#biocLite("ALL")
library(ALL)
#biocLite("hgu95av2.db")
#library(package = affyLib, character.only = TRUE)
#biocLite("Rgraphviz")
library("Rgraphviz")

#tutorial 
# data("geneList")
# data("ALL")
# affyLib <- paste(annotation(ALL), "db", sep = ".")
# library(package = affyLib, character.only = TRUE)


#Load gene GOID data in the format: 
#gene_name    GO:0016021, GO:0016020
geneID2GO <- readMappings("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/gene2go.txt", sep = "\t", IDsep = ",")  

#The 'gene universe' can be all the genes that contain GO annotations in your input 'annotations.txt' file
geneUniverse <- names(geneID2GO) 

# list of interesting genes from a file, just containing a list of interesting genes, one per line
genesOfInterest <- read.table("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/interestinggenes.txt",header=TRUE)
genesOfInterest <- genesOfInterest[which(genesOfInterest$padj < 0.05), ] #all of them are <0.05
genesOfInterest <- as.character(genesOfInterest$gene_name) 

#tell TopGO where these interesting genes appear in the 'geneUniverse' vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
#The geneList object tells TopGO which genes in the gene universe are your genes of interest.

############# Putting the data together into an R object ############
#Here I am using ontology = BP 
myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
#Building most specific GOs .....
#( 77 GO terms found. ) 79
#Build GO DAG topology ..........
# ( 433 GO terms and 867 relations. )444
# Annotating nodes ...............
# ( 90 genes annotated to the GO terms. )97
myMFGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
# Building most specific GOs .....
# ( 85 GO terms found. )88
# 
# Build GO DAG topology ..........
# ( 226 GO terms and 290 relations. )230
# 
# Annotating nodes ...............
# ( 95 genes annotated to the GO terms. )102

#Look at the object 
myGOdata
# ------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  My project 
# 
# Ontology:
#   -  BP 
# 
# 139 available genes (all genes from the array):146
#   - symbol:  MSTRG.15820 MSTRG.9319 MSTRG.9372 MSTRG.931 MSTRG.9378  ...
# - 43  significant genes. 50
# 
# 90 feasible genes (genes that can be used in the analysis):97
#   - symbol:  MSTRG.15820 MSTRG.9319 MSTRG.9372 MSTRG.931 MSTRG.9378  ...
# - 22  significant genes. 29
# 
# GO graph (nodes with at least  1  genes):
#   - a graph with directed edges
# - number of nodes = 433 
# - number of edges = 867 
# 
# ------------------------- topGOdata object -------------------------
myMFGOdata
# ------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  My project 
# 
# Ontology:
#   -  MF 
# 
# 139 available genes (all genes from the array):146
#   - symbol:  MSTRG.15820 MSTRG.9319 MSTRG.9372 MSTRG.931 MSTRG.9378  ...
# - 43  significant genes. 50
# 
# 95 feasible genes (genes that can be used in the analysis):102
#   - symbol:  MSTRG.15820 MSTRG.939 MSTRG.9373 MSTRG.930 MSTRG.28296  ...
# - 29  significant genes. 36
# 
# GO graph (nodes with at least  1  genes):
#   - a graph with directed edges
# - number of nodes = 226 230
# - number of edges = 290 295
# 
# ------------------------- topGOdata object -------------------------
#Ontology = CC
myCCGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
# Building most specific GOs .....
# ( 31 GO terms found. )
# 
# Build GO DAG topology ..........
# ( 137 GO terms and 290 relations. )
# 
# Annotating nodes ...............
# ( 85 genes annotated to the GO terms. )92
myCCGOdata
# ------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  My project 
# 
# Ontology:
#   -  CC 
# 
# 139 available genes (all genes from the array):146
#   - symbol:  MSTRG.15820 MSTRG.9319 MSTRG.9372 MSTRG.931 MSTRG.9378  ...
# - 43  significant genes. 50
# 
# 85 feasible genes (genes that can be used in the analysis):92
#   - symbol:  MSTRG.9319 MSTRG.931 MSTRG.9378 MSTRG.930 MSTRG.23191  ...
# - 27  significant genes. 34
# 
# GO graph (nodes with at least  1  genes):
#   - a graph with directed edges
# - number of nodes = 137 
# - number of edges = 290 
# 
# ------------------------- topGOdata object -------------------------

#Access significant genes
#BP
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata) #22 sig genes 29
#MF
sg_mf <- sigGenes(myMFGOdata)
str(sg_mf)
numSigGenes(myMFGOdata) #29 sig genes 36
#CC
sg_cc <- sigGenes(myCCGOdata)
str(sg_cc)
numSigGenes(myCCGOdata) #27 sig genes 34

############ Performing enrichment tests ############ 
#The p-values computed by the runTest function are unadjusted for multiple testing.
#Fishers exact test
#ontology = BP
resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher") 
#Selecting 'algorithm=classic' means that the GO hierarchy isn't taken into account, 
#so each GO term is tested independently.
resultFisher #The p-values have not been corrected for multiple testing.
# Description: My project 
# Ontology: BP 
# 'classic' algorithm with the 'fisher' test
# 433 GO terms scored: 1 terms with p < 0.01 444
# Annotation data:
#   Annotated genes: 90 97
# Significant genes: 22 29
# Min. no. of genes annotated to a GO: 1 
# Nontrivial nodes: 188 205
#ontology = MF
resultFisherMF <- runTest(myMFGOdata, algorithm="classic", statistic="fisher") 
resultFisherMF #The p-values have not been corrected for multiple testing.
# Description: My project 
# Ontology: MF 
# 'classic' algorithm with the 'fisher' test
# 226 GO terms scored: 0 terms with p < 0.01 230
# Annotation data:
#   Annotated genes: 95 102
# Significant genes: 29 36
# Min. no. of genes annotated to a GO: 1 
# Nontrivial nodes: 121 128
#ontology = CC
resultFisherCC <- runTest(myCCGOdata, algorithm="classic", statistic="fisher") 
resultFisherCC
# Description: My project 
# Ontology: CC 
# 'classic' algorithm with the 'fisher' test
# 137 GO terms scored: 0 terms with p < 0.01 137
# Annotation data:
#   Annotated genes: 85 92
# Significant genes: 27 34
# Min. no. of genes annotated to a GO: 1 
# Nontrivial nodes: 86 88

#Kolmogorov-Smirnov test
#resultKS <- runTest(myGOdata, algorithm = "classic", statistic = "ks") #classic
#resultKS
# Description: My project 
# Ontology: BP 
# 'classic' algorithm with the 'ks' test
# 433 GO terms scored: 0 terms with p < 0.01
# Annotation data:
#   Annotated genes: 90 
# Significant genes: 22 
# Min. no. of genes annotated to a GO: 1 
# Nontrivial nodes: 433 
#resultKS.elim <- runTest(myGOdata, algorithm = "elim", statistic = "ks") #elim
#resultKS.elim
# Description: My project 
# Ontology: BP 
# 'elim' algorithm with the 'ks : 0.01' test
# 433 GO terms scored: 0 terms with p < 0.01
# Annotation data:
#   Annotated genes: 90 
# Significant genes: 22 
# Min. no. of genes annotated to a GO: 1 
# Nontrivial nodes: 433 

#list the top 22 significant results found
#The GenTable function returns a data frame containing the top topNodes GO terms identified by the elim algorithm,
#allRes <- GenTable(myGOdata, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 22)

#BP
FisherRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 22)
write.csv(FisherRes, file = "FisherResBPall.csv")
#MF
FisherResMF <- GenTable(myMFGOdata, classicFisher = resultFisherMF, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 22)
write.csv(FisherResMF, file = "FisherResMFall.csv")
#CC
FisherResCC <- GenTable(myCCGOdata, classicFisher = resultFisherCC, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 22)
write.csv(FisherResCC, file = "FisherResCCall.csv")

#showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
#printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

#The number of GO terms in the TopGO subset of the GO hierarchy (the GO terms annotated to genes in the gene 'universe' input file, 
#plus the ancestors of those GO terms in the GO hierarchy) can be found using:
length(usedGO(myGOdata)) #433 : agrees with the number of nodes in the myGOdata obj. 

pValue.classic <- score(resultFisher)
pValue.classic
pValue.KS <- score(resultKS)[names(pValue.classic)]
gstat <- termStat(myGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
#Found the colMap function here:https://github.com/davetgerrard/LiverProteins/blob/master/scripts/topGO_vignetteTrial.R
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.KS, xlab = "p-value classic", ylab = "p-value KS",main = "p-value differences between Fisher and KS",pch = 19, cex = gSize, col = gCol)
#Some GO terms found significant by the Fisher method are less significant in the KS, as expected.
#The size of the dot is proportional with the number of annotated genes for the respective GO term and its coloring represents the number of significantly 
#differentially expressed genes, with the dark red points having more genes then the yellow ones.

#Lets see how the significant GO terms are distributed over the GO graph.
#subgraph induced by the 5 most significant GO terms as identified by the each algorithm. 
#Significant nodes are represented as rectangles.
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
showSigOfNodes(myGOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')
showSigOfNodes(myGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

