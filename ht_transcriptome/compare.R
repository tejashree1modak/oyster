## get common gene names between two files ##
file1 <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design5/Gene_dfSig_des5_con_trans_12vs5_0.1.csv", 
                  sep = ",", header = TRUE)
genes1 <- file1[,1]
file2 <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design5/Gene_dfSig_des5_con_trans_12vs16_0.1.csv", 
                  sep = ",", header = TRUE)
genes2 <- file2[,1]
common <- genes2[genes2 %in% genes1]
head(common)
common
file1[file1 %in% file2[,1]]
intersect(genes1,genes2)
intersect(file1[,1], file2[,1])
########################################################################################################

### Get annotations for comparisons between treatments #####
library(dplyr)
s4ri <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/deseq/trans_1/common_s4ri.txt", 
                header = F)
s4re <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/deseq/trans_1/common_s4re.txt", 
                 header = F)
s4reri <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/deseq/trans_1/common_s4reri.txt", 
                   header = F)
s4vscon_annot <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/blast/trans_1/S4vsCon_deg_sorted2.csv", 
                          sep = ",", header = TRUE)
s4ri_exclusive <- anti_join(s4ri, s4reri)

colnames(s4re_exclusive) <- c('Iteration_query.def')
colnames(s4ri_exclusive) <- c('Iteration_query.def')
colnames(s4reri) <- c('Iteration_query.def')
colnames(s4re) <- c('Iteration_query.def')
s4re_exclusive <- anti_join(s4re, s4reri)
s4ri_exclusive_annot <- left_join(s4ri_exclusive, s4vscon_annot, by = 'Iteration_query.def')
head(s4ri_exclusive_annot)
write.csv(s4ri_exclusive_annot, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/blast/trans_1/s4ri_exclusive_annot.csv",
          row.names = F)
s4reri_annot <- left_join(s4reri, s4vscon_annot, by = 'Iteration_query.def')
write.csv(s4reri_annot, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/blast/trans_1/s4reri_annot.csv",
          row.names = F)
s4re_exclusive_annot <- left_join(s4re_exclusive, s4vscon_annot, by = 'Iteration_query.def')
write.csv(s4re_exclusive_annot, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/blast/trans_1/s4re_exclusive_annot.csv",
          row.names = F)

########################################################################################################

## map go ids and names to degs ## 
revscon_go_annot <- read.csv("/Users/tejashree/b2gWorkspace/trans2017_final/trans_1/REvsCon/revscon_go_annot.txt", 
                             sep = "\t", header = TRUE)
revscon_deg_sorted <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/blast/trans_1/REvsCon_deg_sorted.csv", 
                               sep = ",", header = TRUE)
revscon_deg_go <- left_join(revscon_deg_sorted, revscon_go_annot, by = 'Iteration_query.def')
write.csv(revscon_deg_go, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/blast/trans_1/revscon_deg_go.csv",
          row.names = F)

