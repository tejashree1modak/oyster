#This script is used to get the top 10 blast hit per transcript to get annotation for genes with 'uncharacterized' annotation
library(plyr)

setwd("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/")

#The 2 csv files below are made on bluewaves by using the python script blast_extract_anno.py.
#Instructions to generate them are in the log section "Uncharacterized".
testx <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/trans_1/S4vsCon_combined_annot_top10.csv", 
                  sep=",", header=TRUE)
testy <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/trans_1/S4vsCon_combined_unchar.csv", 
                   sep=",", header=TRUE)
testnew <- match_df(testx, testy, on = "Iteration_query.def")

#Get the length and unique count of ids
sapply(testnew, function(x) length(unique(x)))
#Remove hits that contains uncharacterized/hypothetical in their annotation
testnew_notunchar <- testnew[!grepl("uncharacterized", testnew$Hit_def),]
sapply(testnew_notunchar, function(x) length(unique(x)))
testnew_notunchar <- testnew_notunchar[!grepl("Uncharacterized", testnew_notunchar$Hit_def),]
sapply(testnew_notunchar, function(x) length(unique(x)))
testnew_notunchar <- testnew_notunchar[!grepl("hypothetical", testnew_notunchar$Hit_def),]
sapply(testnew_notunchar, function(x) length(unique(x)))
#Keep top hit_def and remove others per id
testnew_notunchar_uniq <- testnew_notunchar[!duplicated(testnew_notunchar$Iteration_query.def),]
#making sure we keep all hits and not
sapply(testnew_notunchar_uniq, function(x) length(unique(x)))
#write output 
write.csv(testnew_notunchar_uniq, file = "S4vsCon_trans1_newannot.csv",
          row.names = F)

#Combine new annot for all treatments#
new_annot1 <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/trans_1/REvsCon_trans1_newannot.csv", 
                        sep=",", header=TRUE)
new_annot2 <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/trans_1/RIvsCon_trans1_newannot.csv", 
                        sep=",", header=TRUE)
new_annot3 <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/trans_1/S4vsCon_trans1_newannot.csv", 
                        sep=",", header=TRUE)
new_annot4 <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/trans_2/RI_2vsCon_trans2_newannot.csv", 
                        sep=",", header=TRUE)
new_annot5 <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/trans_2/S4_2vsCon_trans2_newannot.csv", 
                        sep=",", header=TRUE)
new_annot1$treatment <- "revscon_trans1"
new_annot2$treatment <- "rivscon_trans1"
new_annot3$treatment <- "s4vscon_trans1"
new_annot4$treatment <- "ri_2vscon_trans1"
new_annot5$treatment <- "s4_2vscon_trans1"
new_annot_comb <- rbind(new_annot1,new_annot2,new_annot3,new_annot4,new_annot5)
write.csv(new_annot_comb, file = "trans1_2_new_annot_comb.csv",
          row.names = F)