
library(tidyverse)
library(pheatmap)
library(limma)
library(edgeR)

df <- read.csv('gene_exp_add_all_diff.csv',sep = ',',header = T,check.names = F)
colnames(df)
#deplete NA
df <- na.omit(df)

df <- df[!duplicated(df$gene_name),]
rownames(df)<- df$gene_name

df_count <-df[,grepl('count',colnames(df))]
df_FPKM <-df[,grepl('FPKM',colnames(df))]

sum(rowSums(df_count) == 0)
df_count <- df_count[rowSums(df_count) != 0,]
sum(is.na(df_count))

rownames(df_count) <- gsub('-','_',rownames(df_count))

# Find all column names starting with a number and add 'gene' before them
rownames(df_count)[grepl("^\\d", rownames(df_count))] <- paste("gene", rownames(df_count)[grepl("^\\d", rownames(df_count))] ,sep='_')


sum(rowSums(df_FPKM) == 0)
df_FPKM <- df_FPKM[rowSums(df_FPKM) != 0,]
sum(is.na(df_FPKM))

rownames(df_FPKM) <- gsub('-','_',rownames(df_FPKM))

rownames(df_FPKM)[grepl("^\\d", rownames(df_FPKM))] <- paste("gene", rownames(df_FPKM)[grepl("^\\d", rownames(df_FPKM))] ,sep='_')



df_diff_T1 <- read.csv('1.diff_g1g0_T1.csv',header = T,check.names = F,row.names = 1)
df_diff_T1 <- df_diff_T1[abs(df_diff_T1$logFC) >= 2 & df_diff_T1$adj.P.Val < 0.05,]

df_diff_T2 <- read.csv('2.diff_g1g0_T2.csv',header = T,check.names = F,row.names = 1)
df_diff_T2 <- df_diff_T2[abs(df_diff_T2$logFC) >= 2 & df_diff_T2$adj.P.Val < 0.05,]

df_diff_T3 <- read.csv('3.diff_g1g0_T3.csv',header = T,check.names = F,row.names = 1)
df_diff_T3 <- df_diff_T3[abs(df_diff_T3$logFC) >= 2 & df_diff_T3$adj.P.Val < 0.05,]

df_g <- read.csv('CM310.csv',header = T)
df_g$T1_count <- paste0(df_g$pre,'_count')
df_g$T2_count <- paste0(df_g$Day14,'_count')
df_g$T3_count <- paste0(df_g$Day28,'_count')
df_g$T1_FPKM  <- paste0(df_g$pre,'_FPKM')
df_g$T2_FPKM <- paste0(df_g$Day14,'_FPKM')
df_g$T3_FPKM <- paste0(df_g$Day28,'_FPKM')

sample_T2_01 <- c(df_g[df_g$group==0,'T2_count'],df_g[df_g$group==1,'T2_count'])[c(df_g[df_g$group==0,'T2_count'],df_g[df_g$group==1,'T2_count']) %in% colnames(df_count)]
FPKM_T2_01 <- c(df_g[df_g$group==0,'T2_FPKM'],df_g[df_g$group==1,'T2_FPKM'])[c(df_g[df_g$group==0,'T2_FPKM'],df_g[df_g$group==1,'T2_FPKM']) %in% colnames(df_FPKM)]

df_count_T2 <- df_count[rownames(df_diff_T2),sample_T2_01]
df_FPKM_T2 <- df_FPKM[rownames(df_diff_T2),FPKM_T2_01]

df_count_T3 <- df_count[rownames(df_diff_T3),colnames(df_count)[grepl('-03_',colnames(df_count))]]

paste0(df_g$Day14,'_count') %in% colnames(df_count_T2)
colnames(df_count_T2) %in% paste0(df_g$Day14,'_count')
colnames(df_count_T3) %in% paste0(df_g$Day28,'_count')
colnames(df_count)
diff_gene_T2T3 <- intersect(rownames(df_diff_T2),rownames(df_diff_T3))
diff_gene_T2 <- rownames(df_diff_T2)[! rownames(df_diff_T2) %in% diff_gene_T2T3]
diff_gene_T3 <- rownames(df_diff_T3)[! rownames(df_diff_T3) %in% diff_gene_T2T3]
df_count_T2T3 <- df_count[c(diff_gene_T2,diff_gene_T2T3,diff_gene_T3),c(colnames(df_count)[grepl('-02_',colnames(df_count))],colnames(df_count)[grepl('-03_',colnames(df_count))])]
c(colnames(df_count)[grepl('-02_',colnames(df_count))],colnames(df_count)[grepl('-03_',colnames(df_count))])

### 
sample_0_T1 <- df_g[df_g$group==0,'T1_count'][ df_g[df_g$group==0,'T1_count']%in% colnames(df_count)]
sample_0_T2 <- df_g[df_g$group==0,'T2_count'][ df_g[df_g$group==0,'T2_count']%in% colnames(df_count)]
sample_0_T3 <- df_g[df_g$group==0,'T3_count'][ df_g[df_g$group==0,'T3_count']%in% colnames(df_count)]
sample_1_T1 <- df_g[df_g$group==1,'T1_count'][ df_g[df_g$group==1,'T1_count']%in% colnames(df_count)]
sample_1_T2 <- df_g[df_g$group==1,'T2_count'][ df_g[df_g$group==1,'T2_count']%in% colnames(df_count)]
sample_1_T3 <- df_g[df_g$group==1,'T3_count'][ df_g[df_g$group==1,'T3_count']%in% colnames(df_count)]
#53,52,53,47,46,47
cumsum(c(53,52,53,47,46,47))
df_count_in <- df_count[c(diff_gene_T2,diff_gene_T2T3,diff_gene_T3),
                        c(sample_0_T1,sample_0_T2,sample_0_T3,
                          sample_1_T1,sample_1_T2,sample_1_T3)]
### 
sample_0_T1 <- df_g[df_g$group==0,'T1_FPKM'][ df_g[df_g$group==0,'T1_FPKM']%in% colnames(df_FPKM)]
sample_0_T2 <- df_g[df_g$group==0,'T2_FPKM'][ df_g[df_g$group==0,'T2_FPKM']%in% colnames(df_FPKM)]
sample_0_T3 <- df_g[df_g$group==0,'T3_FPKM'][ df_g[df_g$group==0,'T3_FPKM']%in% colnames(df_FPKM)]
sample_1_T1 <- df_g[df_g$group==1,'T1_FPKM'][ df_g[df_g$group==1,'T1_FPKM']%in% colnames(df_FPKM)]
sample_1_T2 <- df_g[df_g$group==1,'T2_FPKM'][ df_g[df_g$group==1,'T2_FPKM']%in% colnames(df_FPKM)]
sample_1_T3 <- df_g[df_g$group==1,'T3_FPKM'][ df_g[df_g$group==1,'T3_FPKM']%in% colnames(df_FPKM)]
#53,52,53,47,46,47
df_FPKM_in <- df_FPKM[c(diff_gene_T2,diff_gene_T2T3,diff_gene_T3),
                        c(sample_0_T1,sample_0_T2,sample_0_T3,
                          sample_1_T1,sample_1_T2,sample_1_T3)]

max(df_FPKM_in)
min(df_FPKM_in)
#plot
##
re_heatmap <- df_FPKM_in
re_heatmap[re_heatmap == 0] <- 1
re_heatmap=log2(re_heatmap)

#scale

scale_re <- t(scale(t(re_heatmap)))
max(scale_re )
min(scale_re )
scale_re[scale_re< -2] <- -2
scale_re[scale_re> 2] <- 2
#
normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}
nor_re <- normalization(scale_re)


library(pheatmap)
library(RColorBrewer)
annotation_col <- data.frame(
                             time=c(rep('T1',53),rep('T2',52),rep('T3',53),rep('T1',47),rep('T2',46),rep('T3',47)),
                             group=c(rep(0,158),rep(1,140)),
                             row.names = colnames(nor_re))
ann_colors = list(
  group = c('0'="#969696", '1'="#273992"),
  time = c(T1 = "#80807f",
           T2 = "#327eba",
           T3 = "#e41e25")
  # time = c(T1 = "#bfc0c0",
  #          #T2 = "#327eba",
  #          T2 = "#3f646f",
  #          T3 = "#f7931e")
)

# A4=10
# pdf('pheatmap.new.pdf',width=1.414*A4, height=A4)
pdf('pheatmap.new.pdf',width=30, height=8)
pheatmap(
  #re_heatmap,
  #scale_re,
  nor_re,
  cluster_rows = T,
  cluster_cols = F,
  #gaps_row = c(34,49),
  gaps_col = c(53,105,158,205,251),
  #color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(100),
  color = colorRampPalette(rev(c(brewer.pal(n = 9, name ="RdBu"),brewer.pal(n = 9, name ="YlGnBu")[8:9])))(100),
  #color = colorRampPalette(rev(c('#f7931e','white','#3f646f')))(100),
  fontsize_row  = 5,
  fontsize_col = 3,
  #annotation_row =annotation_row,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  #annotation_names_row = T
  border_color = F,
  cellwidth = 5,
  cellheight =5
)
dev.off()
