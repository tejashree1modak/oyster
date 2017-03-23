## This script reads a csv formatted as Columns:variables and rows observations
#makes a dotplot. Input file needs to be csv or txt file formatted as:
## treatment expression group
## s4        2.05       s4_8days
## ri        1.5        ri_8days
## con       0.5        con_8days
##install and load libraries 
#install.packages("devtools")
library(devtools)
#install_github("easyGgplot2", "kassambara")
library(ggplot2)
library(easyGgplot2)
#ht2 anova:1-15 rows comparison to con_8, 16-30 comparison to con_11, 31-45 comparison to con_15
#ht2 plot:1-14 rows comparison to con_8, 15-28 comparison to con_11, 29-42 comparison to con_15
# qpcr <- read.table("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/qpcr/calc/qpcr_ht2_dotplot.csv",
#                   header = TRUE, stringsAsFactors= FALSE,sep=",")
#ht1anova:1-19 comp with con6_10, 20-37 comp with con6_13, 38-55 comp with con6_17,56-74 comp with con21
#ht1:1-18 comp with con6_10,day5, 19-35 comp with con6_13,day8, 36-52 comp with con6_17day12,53-70 comp with con21day16
qpcr <- read.table("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/qpcr/calc/qpcr_ht1_dotplot.csv",
                   header = TRUE, stringsAsFactors= FALSE,sep=",")
qpcr_frame <- data.frame(qpcr)
ggplot2.dotplot(data=qpcr_frame[1:18,],xName = 'Treatments', yName='Expression',
                groupName = 'Key',addBoxplot = TRUE,dotsize=1.2,
                xShowTickLabel=TRUE,yShowTickLabel=TRUE,
                xTickLabelFont=c(20,"bold","black"),yTickLabelFont=c(20,"bold","black"),
                xtitle="Treatments",ytitle="Expression",
                xShowTitle=TRUE,yShowTitle=TRUE,
                xtitleFont=c(20,"bold","black"),ytitleFont=c(20,"bold","black"),
                mainTitle = 'Comparison of gene expression with Control on Day 5',
                mainTitleFont=c(24,"bold","black"),
                legendBackground=c("darkolivegreen3", 0.5, "solid", "darkolivegreen"),
                legendTextFont=c(20, "bold","black"),
                backgroundColor="darkolivegreen3",removePanelGrid=TRUE,
                binwidth=14)
#summary(qpcr_frame[29:42,])
# results = aov(Expression ~ Treatments, data=qpcr_frame[56:74,])
# TukeyHSD(results, conf.level = 0.95)
# summary(results)
