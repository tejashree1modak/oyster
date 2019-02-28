
library(plyr)
library(tidyr)
library(ggplot2)

setwd("/Users/tejashree/b2gWorkspace/trans2017_final/trans_1/REvsCon/")


#goslim_annot <- read.csv("/Users/tejashree/b2gWorkspace/trans2017_final/trans_1/REvsCon/revscon_goslim_annot.txt", 
#                     sep = "\t", header = TRUE, stringsAsFactors = F)
# goslim_comb <- paste(goslim_annot[,1], collapse = ";")
# goslim_uniq <- unique(paste(goslim_annot[,1], collapse = ";"))

ga <- read.csv2("/Users/tejashree/b2gWorkspace/trans2017_final/trans_1/REvsCon/revscon_goslim_annot.txt", 
                     sep = "\t", skip = 1, col.names = c("GO_Names"),
                  na.strings = c("", "NA")) %>% na.omit %>% 
        mutate(GO_Names = strsplit(as.character(GO_Names), ";")) %>%
        unnest(GO_Names)
ga_c.freq <- table(ga[grep("C:", ga$GO_Names),])
ga_p.freq <- table(ga[grep("P:", ga$GO_Names),])
ga_f.freq <- table(ga[grep("F:", ga$GO_Names),])
barplot(ga_c.freq)
#bp <- barplot(ga_c.freq, axes=FALSE, axisnames = FALSE)
#text(bp, par("usr")[3], labels = ga[grep("C:", ga$GO_Names),], srt = 45, adj = c(1.1, 1.1), xpd = TRUE, cex = .9)
#axis(2)


