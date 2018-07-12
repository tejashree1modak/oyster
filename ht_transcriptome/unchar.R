#This script is used to get the top 10 blast hit per transcript to get annotation for genes with 'uncharacterized' annotation
library(plyr)

#The 2 csv files below are made on bluewaves by using the python script blast_extract_anno.py.
#Instructions to generate them are in the log section "Uncharacterized".
testx <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/trans_2/S4_2vsCon_trans2_combined_annot_top10.csv", 
                  sep=",", header=TRUE)
testy <- read.csv2("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/trans_2/S4_2vsCon_trans2_combined_uncharac.csv", 
                   sep=",", header=TRUE)
testnew <- match_df(testx, testy, on = "Iteration_query.def")

#Get the length and unique count of ids
sapply(testnew, function(x) length(unique(x)))
#Remove hits that contains uncharacterized/hypothetical
testnew_notunchar <- testnew[!grepl("uncharacterized", testnew$Hit_def),]
sapply(testnew_notunchar, function(x) length(unique(x)))
testnew_notunchar <- testnew_notunchar[!grepl("hypothetical", testnew_notunchar$Hit_def),]
sapply(testnew_notunchar, function(x) length(unique(x)))
#Keep top hit_def and remove others per id
testnew_notunchar_uniq <- testnew_notunchar[!duplicated(testnew_notunchar$Iteration_query.def),]
#making sure we keep all hits and not
sapply(testnew_notunchar_uniq, function(x) length(unique(x)))
#write output 
write.csv(testnew_notunchar_uniq, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/uncharacterized/trans_2/S4_2vsCon_trans2_newannot.csv",
          row.names = F)
