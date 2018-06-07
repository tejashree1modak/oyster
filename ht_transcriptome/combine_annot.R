
library(dplyr)
# source("https://bioconductor.org/biocLite.R")
# biocLite("KEGGREST")
# library(KEGGREST)
# head(keggLink("pathway", "crg"))
# head(keggLink("ko", "map00010"))
# tried above but not very useful.


setwd("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations")

### STEP1:
### Combine BLAST and GO annotation and Fold change ###
## map go ids and names to degs ##

## Input: go_annot.txt :This file is the exported output annotations table from BLAST2GO stored in bg2workspace ##
# This file has: Iteration_query.df and sevral other fields like evalue etc. we want GO IDs and GO Names from this table.
go_annot <- read.csv("/Users/tejashree/b2gWorkspace/trans2017_final/trans_1/REvsCon/revscon_go_annot.txt", 
                             sep = "\t", header = TRUE)

## Input: deg_sorted :This file is made by combining DESeq Gene_sigdf.csv file and blast annotations file combined_annot.csv generated on bluewaves 
## Hot to make this file is described in the 'Get BLAST annotations and fold change' section of BLAST in the project log 2017_transcriptome_analysis
# This file has: Iteration_query.df, foldchange, description 
deg_sorted <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/blast/trans_1/REvsCon_deg_sorted.csv", 
                               sep = ",", header = TRUE)
#Make sure headers match, combine using join 
deg_go <- left_join(deg_sorted, go_annot, by = 'Iteration_query.def')

#write to csv if you want but change name or will be overwritten next time #
write.csv(deg_go, file = "/trans_1/deg_go.csv",
          row.names = F)
# At the end of this section we have: Iteration_query.def, fold change, description, GO IDs, GO Names #


####STEP2:
### Get KEGG files ready ###
## Making these input files are described in 'KO to kegg pathway mapping:' section of KEGG annotaions section of the project log. 

## Input: File with pathway num to KO ID mapping ##
pathway2ko <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/kegg/pathway2ko.txt", 
                        sep="\t", header=TRUE) 
## Input: file that has pathway num and name ## 
pathwaynum <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/kegg/pathway_num.txt", 
                        sep="\t", header=TRUE) 
pathwayname_num_ko <- left_join(pathway2ko, pathwaynum, by = 'pathway')
head(pathwayname_num_ko)
ko_pathway_pathname <- pathwayname_num_ko %>% group_by(KO) %>% 
                        summarize(pathway = paste(pathway, collapse = ","), pathway_name = paste(pathway_name, collapse = ","))
head(ko_pathway_pathname)
# At the end of this section we have: KO, pathwaynum, pathwayname #


####STEP3:
## Read in KO annotations output by KAAS ##
## Input: upreg file ##
# output of submission to KAAS annotation server with Iteration_query.def, KO #
ko_up <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/kegg/trans_1/REvsCon_upreg.ko.txt", 
                           sep="\t", header=FALSE)
##Input: downreg file ##
ko_down <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/kegg/trans_1/REvsCon_downreg.ko.txt", 
                           sep="\t", header=FALSE)
# write headers #
colnames(ko_up) <- c('Iteration_query.def','KO')
colnames(ko_down) <- c('Iteration_query.def','KO')
head(ko_up)
head(ko_down)
# remove first row that was added for KAAS:' #revscon_upreg '
ko_up2 <- ko_up[-1,]
ko_down2 <- ko_down[-1,]
# verify that line is gone and store in new variable #
head(ko_up2)
head(ko_down2)
# merge dfs up and down #
ko <- rbind(ko_up2, ko_down2)
# At the end of this section we have: Iteration_query.def, KO #


####STEP4:
## Combine to get Iteration_query.def, KO, pathwaynum, pathwayname ## 
deg_ko <- left_join(ko, ko_pathway_pathname, by = 'KO')
#wite output as csv if wanted #
#Change name if yuo want to save to avoid overwrite!##
write.table(deg_ko, file = "/trans_1/deg_ko.csv",
            row.names = F, quote = FALSE, sep = '\t')
##open in shell and use awk to change gene names to short gene names output as new file##  
##awk -F'\t' -v OFS='\t' '{ split($1, a, /_/); print a[2], $2, $3, $4 }' deg_ko.csv > deg_ko_r.csv ##
#open in vim and add header to Iteration_query.def ##
# At the end of this we have: Iteration_query.def, KO, pathwaynum, pathwayname #


####STEP5:
## Combine all annotations: Iteration_query.def, Foldchange, Description, GO IDs, GO Names, KO, Pathway, pathwayname
## Input: deg_go file :Output of section1- Combine BLAST and GO annotation and Fold change
deg_go <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_1/deg_go.csv", 
                            sep=",", header=TRUE) 
## Input: Output file of section3: Combine to get Iteration_query.def, KO, pathwaynum, pathwayname edited using awk and header aded ##
deg_ko_r <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_1/deg_ko_r.csv", 
                            sep="\t", header=TRUE) 
## Join both ## 
deg_ko_go <- left_join(deg_ko_r, deg_go, by = 'Iteration_query.def')
deg_ko_go <- deg_ko_go[c("Iteration_query.def", "log2FoldChange", "Hit_def", "GO.Names", "pathway_name", "pathway","GO.IDs", "KO" )]

## Write output as csv ## 
#Change name to avoid overwrite!! ##
write.table(deg_ko_go, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_1/deg_ko_go.csv",
            row.names = F, quote = FALSE, sep = '\t')
