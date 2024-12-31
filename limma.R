
library(nlme)
library(emmeans)
setwd('/data/home/teacher/project/tongren_sc/seqrna')
df_hhh <- read.csv('gene_exp_add_all_diff.xls',
                   sep='\t',check.names = F)

colnames(df_hhh)
#deplete NA
df_hhh <- na.omit(df_hhh)
#gene names as rownames
df_hhh <- df_hhh[!duplicated(df_hhh$gene_name),]
rownames(df_hhh)<- df_hhh$gene_name
df_hhh_count <-df_hhh[,grepl('count',colnames(df_hhh))]
rm(df_hhh)
colnames(df_hhh_count)
df_hhh_count <-df_hhh_count[,!grepl('-04_',colnames(df_hhh_count))]
#Remove sum rows with 0 values
sum(rowSums(df_hhh_count) == 0)
df_hhh_count <- df_hhh_count[rowSums(df_hhh_count) != 0,]
sum(is.na(df_hhh_count))

rownames(df_hhh_count) <- gsub('-','_',rownames(df_hhh_count))


library(genefilter)
df_group_310 <- read.table('CM310.csv',sep = ',',header = T,colClasses = "character")
df_group_310$pre <- paste(df_group_310$pre,'count',sep = '_')
df_group_310$Day14 <- paste(df_group_310$Day14,'count',sep = '_')
df_group_310$Day28 <- paste(df_group_310$Day28,'count',sep = '_')
library(reshape2)
df_hhh_g_1 <- reshape2::melt(df_group_310[,c(2:5)],id.vars = 'group',variable.name = 'time_old',value.name = 'sample')

df_hhh_g <-data.frame(sample=colnames(df_hhh_count))
df_hhh_g <- merge(df_hhh_g ,df_hhh_g_1,by='sample')
df_hhh_g$time <- ifelse(grepl('-01_',df_hhh_g$sample),'T1',
                        ifelse(grepl('-02_',df_hhh_g$sample),'T2',
                               ifelse(grepl('-03_',df_hhh_g$sample),'T3','错了')))

rownames(df_hhh_g) <- df_hhh_g$sample
df_hhh_g <- separate(df_hhh_g,col = 'sample',into = c("ID", "new_column2"), sep = "-")
df_hhh_g <- df_hhh_g[,c(1,3,5)]


write.csv(df_hhh_g,'310项目分组信息_new.csv')


##################################310 limma ##################


library(genefilter)

#using for group comparison
df_hhh_g_old = df_hhh_g
#g0,g1 at T1
df_hhh_g=df_hhh_g_old
df_hhh_g$time[df_hhh_g$time=='T1'] = 'T0'
y <- DGEList(counts = df_hhh_count[,rownames(df_hhh_g)])
sum(rownames(df_hhh_g) != colnames(y))
#Keep genes with total counts more than 50.
A <- rowSums(y$counts)
isexpr <- A > 50
y <- y[isexpr, , keep.lib.size = FALSE]
dim(y)
#Apply scale normalization:
y <- calcNormFactors(y)
design <- model.matrix(~  time*group, data =df_hhh_g)
v <- voom(y, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
summary(decideTests(fit))
diff_g1g0_T1 <- topTable(fit, coef = "group1", n = Inf)
write.csv(diff_g1g0_T1,'diff_g1g0_T1.csv')

#g0,g1 at T2
df_hhh_g=df_hhh_g_old
df_hhh_g$time[df_hhh_g$time=='T2'] = 'T0'
y <- DGEList(counts = df_hhh_count[,rownames(df_hhh_g)])
sum(rownames(df_hhh_g) != colnames(y))

#Keep genes with total counts more than 50.
A <- rowSums(y$counts)
isexpr <- A > 50
y <- y[isexpr, , keep.lib.size = FALSE]
dim(y)
#Apply scale normalization:
y <- calcNormFactors(y)
design <- model.matrix(~  time*group, data =df_hhh_g)
v <- voom(y, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
summary(decideTests(fit))
diff_g1g0_T2 <- topTable(fit, coef = "group1", n = Inf)
write.csv(diff_g1g0_T2,'diff_g1g0_T2.csv')

#g0,g1 at T3
df_hhh_g=df_hhh_g_old
df_hhh_g$time[df_hhh_g$time=='T3'] = 'T0'
y <- DGEList(counts = df_hhh_count[,rownames(df_hhh_g)])
sum(rownames(df_hhh_g) != colnames(y))
#Keep genes with total counts more than 50.
A <- rowSums(y$counts)
isexpr <- A > 50
y <- y[isexpr, , keep.lib.size = FALSE]
dim(y)
#Apply scale normalization:
y <- calcNormFactors(y)

design <- model.matrix(~  time*group, data =df_hhh_g)
v <- voom(y, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
summary(decideTests(fit))
###
diff_g1g0_T3 <- topTable(fit, coef = "group1", n = Inf)
write.csv(diff_g1g0_T3,'diff_g1g0_T3.csv')

#T2,T3 at group 1
df_hhh_g=df_hhh_g_old
df_hhh_g$group[df_hhh_g$group=='0'] = '2'
###
y <- DGEList(counts = df_hhh_count[,rownames(df_hhh_g)])
sum(rownames(df_hhh_g) != colnames(y))

#Keep genes with total counts more than 50.
A <- rowSums(y$counts)
isexpr <- A > 50
y <- y[isexpr, , keep.lib.size = FALSE]
dim(y)
#Apply scale normalization:
y <- calcNormFactors(y)

design <- model.matrix(~  time*group, data =df_hhh_g)
v <- voom(y, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
summary(decideTests(fit))
###
diff_T3T1_g1 <- topTable(fit, coef = "timeT3", n = Inf)
diff_T2T1_g1 <- topTable(fit, coef = "timeT2", n = Inf)
write.csv(diff_T3T1_g1,'diff_T3T1_g1.csv')
write.csv(diff_T2T1_g1,'diff_T2T1_g1.csv')

#T2,T3 at group 0
df_hhh_g=df_hhh_g_old
###
y <- DGEList(counts = df_hhh_count[,rownames(df_hhh_g)])
sum(rownames(df_hhh_g) != colnames(y))

#Keep genes with total counts more than 50.
A <- rowSums(y$counts)
isexpr <- A > 50
y <- y[isexpr, , keep.lib.size = FALSE]
dim(y)
#Apply scale normalization:
y <- calcNormFactors(y)

design <- model.matrix(~  time*group, data =df_hhh_g)
v <- voom(y, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
summary(decideTests(fit))
diff_T3T1_g0 <- topTable(fit, coef = "timeT3", n = Inf)
diff_T2T1_g0 <- topTable(fit, coef = "timeT2", n = Inf)
write.csv(diff_T3T1_g0,'diff_T3T1_g0.csv')
write.csv(diff_T2T1_g0,'diff_T2T1_g0.csv')



