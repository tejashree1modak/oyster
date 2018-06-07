library("dplyr")
library(ggplot2)
install.packages("pca3d")
library(pca3d)
install.packages("rgl")
library(rgl)
install.packages("devtools")
library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)
install.packages("ggfortify")
library(ggfortify)
install.packages("VennDiagram")
library(VennDiagram)
library(gridExtra)
install.packages("UpSetR")
library(UpSetR)
library(plyr)
library(reshape2)

#For mstrg: 
gene_sig <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/Gene_dfSig_mstrg_fd_padj.csv", 
                     sep = ",", header = TRUE)
annot <- read.csv("/Users/tejashree/b2gWorkspace/design1/Gene_dfSig_mstrg_annot.csv", sep = ",", header = TRUE)
m <- merge(gene_sig, annot)
head(m)
write.csv(m, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/de_mstrg.csv")

#For geneids:
geneids_gene_sig <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/Gene_dfSig_geneids_fd_padj.csv", 
                             sep = ",", header = TRUE)
geneids_annot <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/Gene_dfSig_geneids_annot.csv", 
                          sep = ",", header = TRUE)
geneids_m <- merge(geneids_gene_sig,geneids_annot)
head(geneids_m)
write.csv(geneids_m, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/de_geneids.csv")


knownupdwn <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/de_mstrg_sorted_known_genes_up_down.csv", 
                  sep = ",", header = TRUE)
##################
#Get annotations + foldchange csv
anno <- read.csv("/tmp/3", 
                 sep = ",", header = TRUE)
fldchange <- read.csv("/tmp/1", 
                      sep = ",", header = TRUE)
m <- merge(fldchange, anno)
head(m)
write.csv(m, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2017/transcriptome_2017/data_analysis/blast/trans_1/S4vsCon_deg_sorted2.csv",
          row.names = F)

###########
#design2
##########
gene_sig_des2 <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design2/Gene_dfSig_fd_padj_des2.csv", 
                     sep = ",", header = TRUE)
annot_des2 <- read.csv("/Users/tejashree/b2gWorkspace/design2/Genedf_sig_annot.csv", sep = ",", header = TRUE)
m <- merge(gene_sig_des2, annot_des2)
head(m)
write.csv(m, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design2/de_des2.csv")

########
#design4
########
gff_annot <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design4/Gene_dfSig_des4_5d_gff_annotation.csv",
                  sep = ",", header = TRUE)
uniq_gff_annot <- gff_annot[!duplicated(gff_annot[,c('gene_id')]),]
head(uniq_gff_annot)
write.csv(uniq_gff_annot, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design4/Gene_dfSig_des4_5d_gff_annotation_uniq.csv")
#merge
gene_sig_des4 <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design4/Gene_dfSig_des4_5d_fd_padj.csv",
                          sep = ",", header = TRUE)
gene_ids_annot_des4 <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design4/Gene_dfSig_des4_5d_gff_annotation_uniq.csv",
                                sep = ",", header = TRUE)
m_des4 <- merge(gene_sig_des4,gene_ids_annot_des4)
head(m_des4)
write.csv(m_des4, file = "/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design4/de_des4.csv",
          row.names=FALSE, quote = FALSE)
###################################
#PCA2 plot for genecounts 
###################################
t <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/de_genecount_trans.csv")
t_sur <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/de_genecount_trans_survival.csv")
#Useful link for visualization of the pca: 
#https://tgmstat.wordpress.com/2013/11/28/computing-and-visualizing-pca-in-r/
#only genes
geneid <- t[,1]
pca <- prcomp(t[,-1], scale.=TRUE)
print(pca)
summary(pca)
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, 
              ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
biplot(pca)
#genes&survival
t_sur_slice <- t_sur[,90:104]
gid <- c('con5d', 'con12d', 'con16d', 'rif5d', 'rif12d', 'rif16d')
gid
t_sur_slice$geneid <- gid
t_sur_slice[,-16]
pca_sur <- prcomp(t_sur[,-1], scale.=TRUE)
pca_sur_slice <- prcomp(t_sur_slice[,-16], scale.=TRUE)
summary(pca_sur)
summary(pca_sur_slice)
print(pca_sur)
print(pca_sur_slice)
biplot(pca_sur)
biplot(pca_sur_slice)
#plot below used the info from:
#https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
autoplot(prcomp(t[,-1], scale.=TRUE), data = t, colour = 'gene_id')

# -----trial for group-----
aa <- read.csv("/tmp/a.csv")
aaa <- as.data.frame(aa)
#not worked yet! 

##MDS plot 
d <- dist(t_sur)
d_t <- dist(t)
fit <- cmdscale(d,eig=TRUE, k=2, x.ret = TRUE) # k is the number of dimensions
fit_t <- cmdscale(d_t, eig=TRUE, k=2, x.ret = TRUE) # k is the number of dimensions
fit # view results
fit_t

# plot solution with survival  
x_fit <- fit$points[,1]
y_fit <- fit$points[,2]
plot(x_fit, y_fit, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS")
labels = c("con5d", "con12d","cpn16d","rif5d", "rif12d", "rif16d")
text(x_fit, y_fit, labels = labels, cex=1, pos = 1)
# plot solution without survival 
x_t_fit <- fit_t$points[,1]
y_t_fit <- fit_t$points[,2]
plot(x_t_fit, y_t_fit, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS")
labels = c("con5d", "con12d","con16d","rif5d", "rif12d", "rif16d")
text(x_t_fit, y_t_fit, labels = labels, cex=1, pos = 1)
# look up pcoa function, you have negative eigenvalues can solve that by using a correction in pcoa function.

########################################################
#Venn Diagram
########################################################
#Venn diagram for design 1
#resource: https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html
v <- read.csv("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/venn.csv", header = T)
con <- select(v, con5d, con12d, con16d)
treat <- select(v, rif5d, rif12d, rif16d)
con_overlap <- calculate.overlap(con = list("Control" = con))
con5d <- nrow(subset(v, con5d == 1))
con12d <- nrow(subset(v, con12d == 1))
con16d <- nrow(subset(v, con16d == 1))
con512 <- nrow(subset(v, con5d == 1 & con12d == 1))
con1216 <- nrow(subset(v, con12d == 1 & con16d == 1))
con516 <- nrow(subset(v, con5d == 1 & con16d == 1))
con51216 <- nrow(subset(v, con5d == 1 & con12d == 1 & con16d == 1))
venna <- draw.triple.venn(area1 = con5d, area2 = con12d, area3 = con16d, n12 = con512, n23 = con1216, n13 = con516, n123 = con51216,
                 category = c("con5d", "con12d", "con16d"), fill = c("skyblue", "pink1", "mediumorchid"))
grid.newpage()
rif5d <- nrow(subset(v, rif5d == 1))
rif5d
rif12d <- nrow(subset(v, rif12d == 1))

rif16d <- nrow(subset(v, rif16d == 1))
rif512 <- nrow(subset(v, rif5d == 1 & rif12d == 1))
rif1216 <- nrow(subset(v, rif12d == 1 & rif16d == 1))
rif516 <- nrow(subset(v, rif5d == 1 & rif16d == 1))
rif51216 <- nrow(subset(v, rif5d == 1 & rif12d == 1 & rif16d == 1))
rif51216
vennb <- draw.triple.venn(area1 = rif5d, area2 = rif12d, area3 = rif16d, n12 = rif512, n23 = rif1216, n13 = rif516, n123 = rif51216,
                 category = c("rif5d", "rif12d", "rif16d"), fill = c("skyblue", "pink1", "mediumorchid"))
grid.arrange(gTree(children=venna), gTree(children=vennb), ncol=2)

grid.newpage()
conrif5d <- nrow(subset(v, con5d ==1 & rif5d == 1))
conrif12d <- nrow(subset(v, con12d ==1 & rif12d == 1))
conrif16d <- nrow(subset(v, con16d ==1 & rif16d == 1))
venn1 <- draw.pairwise.venn(con5d, rif5d, conrif5d, category = c("con5d", "rif5d"), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                   cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = FALSE)
venn2 <- draw.pairwise.venn(con12d, rif12d, conrif12d, category = c("con12d", "rif12d"), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                   cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = FALSE)
venn3 <- draw.pairwise.venn(con16d, rif16d, conrif16d, category = c("con16d", "rif16d"), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                   cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = FALSE)
grid.arrange(gTree(children=venn1), gTree(children=venn2),
             gTree(children=venn3), ncol=3)

##Venn diagram for design1_trans: Using transcript count. 
trans_ct <- read.csv2("~/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/transcript_count_matrix_names.csv",
                       sep = ",", header = TRUE)
degs_trans <- read.csv2("~/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/trans/Gene_dfSig_des1_trans.csv",
                        sep = ",", header = TRUE) 
trans_ct_degs <- subset(trans_ct, transcript_id %in% degs_trans$X) #Can also do with dplyr: trans_ct %>% filter(transcript_id %in% degs_trans$X)
# can de done this way too: trans_ct_degs$con5d = ifelse(trans_ct_degs$con5d >= 1 , 1, 0)
binary <- trans_ct_degs %>% mutate(con5d = ifelse(con5d >= 1, 1, 0), 
                         con12d = ifelse(con12d >= 1, 1, 0),
                         con16d = ifelse(con16d >= 1, 1, 0),
                         rif5d = ifelse(rif5d >= 1, 1, 0),
                         rif12d = ifelse(rif12d >= 1, 1, 0),
                         rif16d = ifelse(rif16d >= 1, 1, 0)) 
write.csv(binary, file = "~/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/trans/venn.csv",
            row.names = F, quote = F)
             
con5d_list <- trans_ct_degs$transcript_id[which(trans_ct_degs$con5d == 1)] 
con12d_list <- trans_ct_degs$transcript_id[which(trans_ct_degs$con12d == 1)] 
con16d_list <- trans_ct_degs$transcript_id[which(trans_ct_degs$con16d == 1)] 
rif5d_list <- trans_ct_degs$transcript_id[which(trans_ct_degs$rif5d == 1)] 
rif12d_list <- trans_ct_degs$transcript_id[which(trans_ct_degs$rif12d == 1)] 
rif16d_list <- trans_ct_degs$transcript_id[which(trans_ct_degs$rif16d == 1)] 

upset(binary, nsets = 6, nintersects = NA)

#Trans2017_final analysis
venn_2017_vsCon <- draw.triple.venn(area1 = 1550, area2 = 1534, area3 = 2269, n12 = 774, n23 = 802, n13 = 787, n123 = 540,
                          category = c("RIvsCon", "REvsCon", "S4vsCon"), fill = c("skyblue", "yellow", "lightgreen"), scaled = TRUE, euler.d = TRUE)
grid.draw(venn_2017_vsCon)
dev.copy(pdf, "Venn_2017_vsCon.pdf")
dev.off()
venn_2017_vsProb <- draw.pairwise.venn(area1 = 1229 , area2 = 1757, cross.area = 505, scaled = TRUE, euler.d = TRUE, fill = c("skyblue", "lightgreen"),
                                       category = c("RIvsRE", "S4vsRE"))
grid.draw(venn_2017_vsProb)
dev.copy(pdf, "Venn_2017_vsProb.pdf")
dev.off()