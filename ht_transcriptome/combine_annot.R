
library(dplyr)
# source("https://bioconductor.org/biocLite.R")
# biocLite("KEGGREST")
# library(KEGGREST)
# head(keggLink("pathway", "crg"))
# head(keggLink("ko", "map00010"))
# tried above but not very useful.


setwd("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/")

### STEP1:
### Combine BLAST and GO annotation and Fold change ###
## map go ids and names to degs ##

## Input: go_annot.txt :This file is the exported output annotations table from BLAST2GO stored in bg2workspace ##
# This file has: Iteration_query.df and sevral other fields like evalue etc. we want GO IDs and GO Names from this table.
go_annot <- read.csv("/Users/tejashree/b2gWorkspace/trans2017_final/trans_2/S4_2vsCon/s4_2vscon_go_annot.txt", 
                             sep = "\t", header = TRUE)

## Input: deg_sorted :This file is made by combining DESeq Gene_sigdf.csv file and blast annotations file combined_annot.csv generated on bluewaves 
## How to make this file is described in the 'Get BLAST annotations and fold change' section of BLAST in the project log 2017_transcriptome_analysis
# This file has: Iteration_query.df, foldchange, description 
deg_sorted <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/blast/trans_2/S4_2vsCon_deg_sorted.csv", 
                               sep = ",", header = TRUE)
#Make sure headers match, combine using join 
deg_go <- left_join(deg_sorted, go_annot, by = 'Iteration_query.def')

#write to csv if you want but change name or will be overwritten next time #
write.csv(deg_go, file = "deg_go.csv",
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
ko_up <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/kegg/trans_2/S4_2vsCon_upreg_ko.txt", 
                           sep="\t", header=FALSE)
##Input: downreg file ##
ko_down <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/kegg/trans_2/S4_2vsCon_downreg_ko.txt", 
                           sep="\t", header=FALSE)
## !! Use this Input code if you dont have separate ko_up/ko_down files !!
# skip next few lines of code in Step3 if you use this since you already have ko dataframe ready
ko <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/kegg/trans_3/", 
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
#Change name if you want to save to avoid overwrite!##
write.table(deg_ko, file = "deg_ko.csv",
            row.names = F, quote = FALSE, sep = '\t')
##open in shell and use awk to change gene names to short gene names output as new file##  
##awk -F'\t' -v OFS='\t' '{ split($1, a, /_/); print a[2], $2, $3, $4 }' deg_ko.csv > deg_ko_r.csv ##
#open in vim and add header to Iteration_query.def ##
# At the end of this we have: Iteration_query.def, KO, pathwaynum, pathwayname #


####STEP5:
## Combine all annotations: Iteration_query.def, Foldchange, Description, GO IDs, GO Names, KO, Pathway, pathwayname
## Input: deg_go file :Output of section1- Combine BLAST and GO annotation and Fold change
deg_go <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_2/deg_go.csv", 
                            sep=",", header=TRUE) 
## Input: Output file of section3: Combine to get Iteration_query.def, KO, pathwaynum, pathwayname edited using awk and header aded ##
deg_ko_r <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_2/deg_ko_r.csv", 
                            sep="\t", header=TRUE) 
## Join both ## 
deg_ko_go <- left_join(deg_ko_r, deg_go, by = 'Iteration_query.def')
deg_ko_go <- deg_ko_go[c("Iteration_query.def", "log2FoldChange", "Hit_def", "GO.Names", "pathway_name", "pathway","GO.IDs", "KO" )]

## Write output as csv ## 
#Change name to avoid overwrite!! ##
write.table(deg_ko_go, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_2/deg_ko_go.tsv",
            row.names = F, quote = FALSE, sep = '\t')
#Change name of output file and can delete intermediate files # 

#STEP6:
#### Add column for treatment with same rows for deg_ko_go.csv ###
#after data for each treatment has gone through above steps #
#this can be used to combine all deg_ko_go tables into one table to see shared genes#
deg_ko_go_1 <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_2/ri_2vscon_deg_ko_go.tsv", 
                      sep="\t", header=TRUE)
deg_ko_go_2 <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_2/s4_2vscon_deg_ko_go.tsv", 
                         sep="\t", header=TRUE)
#deg_ko_go_3 <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_1/rivscon_deg_ko_go.tsv", 
#                         sep="\t", header=TRUE)
deg_ko_go_1$treatment <- "r1_2vscon_trans2"
deg_ko_go_2$treatment <- "s4_2vscon_trans2"
#deg_ko_go_3$treatment <- "rivscon_trans1"
#deg_ko_go_comb <- rbind(deg_ko_go_1,deg_ko_go_2,deg_ko_go_3)
deg_ko_go_comb <- rbind(deg_ko_go_1,deg_ko_go_2)
write.table(deg_ko_go_comb, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_2/deg_ko_go_comb.tsv",
            row.names = F, quote = FALSE, sep = '\t')
#open deg_ko_go.tsv using excel#
#color rows by treatment#
#sort first by pathway_name then GO_names then log2foldchange#
#save as : xlsb in curated dir #

#Step7#
#combining trans1 and trans2 combined annot#
deg_ko_go_trans1 <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_1/reris4vscon_deg_ko_go_comb.tsv", 
                         sep="\t", header=TRUE)
deg_ko_go_trans2 <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans_2/ris4_2vscon_deg_ko_go_comb.tsv", 
                         sep="\t", header=TRUE)
deg_ko_go_trans_comb <- rbind(deg_ko_go_trans1,deg_ko_go_trans2)
write.table(deg_ko_go_trans_comb, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/combined_annotations/trans1_2_deg_ko_go_comb.tsv",
            row.names = F, quote = FALSE, sep = '\t')
#open deg_ko_go.tsv using excel#
#color rows by treatment#
#sort first by pathway_name then GO_names then log2foldchange#
#save as : xlsb in curated dir #

#Step8
#Pulling out data per question asked
#Obj1: RI vs RE
rivsre <- deg_ko_go_trans_comb[deg_ko_go_trans_comb$treatment == 'rivscon_trans1' | deg_ko_go_trans_comb$treatment == 'revscon_trans1', ]
write.table(rivsre, file = "rivsre_deg_ko_go_comb.tsv",
            row.names = F, quote = FALSE, sep = '\t')
ri <- deg_ko_go_trans_comb[deg_ko_go_trans_comb$treatment == 'rivscon_trans1', ]
re <- deg_ko_go_trans_comb[deg_ko_go_trans_comb$treatment == 'revscon_trans1', ]
#the reason why im using hit_def instead of iteration_id.def is becasue multiple ids can map to the same annotation. 
#I want to narrow down the number to look at manually
rionly <- anti_join(ri, re, by = 'Hit_def')
write.table(rionly, file = "rionly_deg_ko_go_comb.tsv",
            row.names = F, quote = FALSE, sep = '\t')
#also tried using id.
#need to think about which one is better to use. 
rionly_id <- anti_join(ri, re, by = 'Iteration_query.def')

rivsre_reonly <- anti_join(re, ri, by = 'Hit_def')
write.table(rivsre_reonly, file = "reonly_rivsre_deg_ko_go_comb.tsv",
            row.names = F, quote = FALSE, sep = '\t')
#Obj1: S4 vs RE
s4 <- deg_ko_go_trans_comb[deg_ko_go_trans_comb$treatment == 's4vscon_trans1', ]
s4only <- anti_join(s4, re, by = 'Hit_def')
write.table(s4only, file = "s4only_deg_ko_go_comb.tsv",
            row.names = F, quote = FALSE, sep = '\t')
#also tried using id.
s4only_id <- anti_join(s4, re, by = 'Iteration_query.def')

s4vsre_reonly <- anti_join(re, s4, by = 'Hit_def')
write.table(s4vsre_reonly, file = "reonly_s4vsre_deg_ko_go_comb.tsv",
            row.names = F, quote = FALSE, sep = '\t')
#getting ids common in rivsre_reonly and s4vsre_reonly
reonly_s4rivscon <- inner_join(rivsre_reonly, s4vsre_reonly, by = 'Hit_def') #321 #is fewer than venn but generated by comparing hit_def
write.table(reonly_s4rivscon, file = "reonly_s4rivscon_byhitdef_deg_ko_go_comb.tsv",
            row.names = F, quote = FALSE, sep = '\t')
#Obj1: RI vs S4 6h,24h