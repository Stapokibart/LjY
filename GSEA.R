library(ggplot2)
library(clusterProfiler)

save_dir = 'GSEA.result.nodiff.new.20241224'
if (!dir.exists(save_dir)) {
  dir.create(save_dir,recursive=T)
}

file_paths <- c('2.diff_g1g0_T2.csv','3.diff_g1g0_T3.csv','5.diff_T2T1_g1.csv','7.diff_T3T1_g1.csv')

# Read each CSV file in a loop for enrichment analysis
for (file_path in file_paths) {
  data <- read.csv(file_path)
  geneList <- data$logFC
  names(geneList) <- as.character(data$X)  # Ensure that the gene ID is character
  geneList <- sort(geneList, decreasing = TRUE)  # descending order
  kegg_gmt_c2_1 <- read.gmt("GSEA/c2.cp.reactome.v2024.1.Hs.symbols.gmt") 
  kegg_gmt_h <- read.gmt("GSEA/h.all.v2024.1.Hs.symbols.gmt")
  kegg_gmt <- rbind(kegg_gmt_c2_1 ,kegg_gmt_h)
  gsea <- GSEA(geneList,TERM2GENE = kegg_gmt) # GSEA analysis
  a <- gsea@result
  if(nrow(a)==0){
    message("gsea.result.null", file_path)
    next
  }else{
  GESA_result_name  <- sub(".csv", "_GSEA_symbols.csv", basename(file_path))
  write.csv(a, file = file.path(save_dir, GESA_result_name))
  #Visualization
  # dotplot
  dotplot(gsea, showCategory = 10)
  GSEA_dotplot_file_pdf <- sub(".csv", "_gsea_dotplot.pdf", basename(file_path))
  ggsave(file.path(save_dir,GSEA_dotplot_file_pdf),width = 10,height = 10,limitsize = FALSE)
  GSEA_dotplot_file_png <- sub(".csv", "_gesa_dotplot.png", basename(file_path))
  ggsave(file.path(save_dir,GSEA_dotplot_file_png),width = 10,height = 10,limitsize = FALSE,
         device = "png")
  }
}

# ploting
#read data
# T2T3_g01
df_T2_g01 <- read.csv('GSEA.result.nodiff.new.20241224/2.diff_g1g0_T2_GSEA_symbols.csv',row.names = 1)
df_T3_g01 <- read.csv('GSEA.result.nodiff.new.20241224/3.diff_g1g0_T3_GSEA_symbols.csv',row.names = 1)

#picked signaling
sx_T23_g10 <- c(
  'REACTOME_INTERLEUKIN_10_SIGNALING',
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
  'REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING',
  'HALLMARK_IL2_STAT5_SIGNALING',
  'HALLMARK_KRAS_SIGNALING_UP',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING',
  'REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING',
  'REACTOME_NUCLEAR_EVENTS_STIMULATED_BY_ALK_SIGNALING_IN_CANCER',
  'REACTOME_INTERLEUKIN_2_FAMILY_SIGNALING',
  'REACTOME_STAT3_NUCLEAR_EVENTS_DOWNSTREAM_OF_ALK_SIGNALING',
  'REACTOME_G_ALPHA_I_SIGNALLING_EVENTS',
  'REACTOME_CONSTITUTIVE_SIGNALING_BY_ABERRANT_PI3K_IN_CANCER',
  'REACTOME_INTERLEUKIN_RECEPTOR_SHC_SIGNALING',
  'REACTOME_SHC1_EVENTS_IN_EGFR_SIGNALING',
  'REACTOME_INTERLEUKIN_3_INTERLEUKIN_5_AND_GM_CSF_SIGNALING',
  'REACTOME_SIGNALING_BY_INTERLEUKINS',
  'REACTOME_INTERLEUKIN_7_SIGNALING',
  'REACTOME_G_ALPHA_S_SIGNALLING_EVENTS',
  'REACTOME_SHC1_EVENTS_IN_ERBB4_SIGNALING',
  'REACTOME_PURINERGIC_SIGNALING_IN_LEISHMANIASIS_INFECTION'
  
)

df_T2_g01 <- df_T2_g01[sx_T23_g10,]
df_T3_g01 <- df_T3_g01[sx_T23_g10,]
df_T2_g01 <- na.omit(df_T2_g01)
df_T3_g01 <- na.omit(df_T3_g01)

#T2
df_T2_g01_in <- df_T2_g01[,c('ID','NES','p.adjust')]
df_T2_g01_in$Regulation <- ifelse(df_T2_g01_in$NES>0,'Up','Down')
df_T2_g01_in$category <- 'T2'


#T3
df_T3_g01_in <- df_T3_g01[,c('ID','NES','p.adjust')]
df_T3_g01_in$Regulation <- ifelse(df_T3_g01_in$NES>0,'Up','Down')
df_T3_g01_in$category <- 'T3'


df_T23 <- rbind(df_T2_g01_in,df_T3_g01_in)
df_T23$p.adjust.log10 <- -log10(df_T23$p.adjust)
df_T23$ID <- tolower(df_T23$ID)

library(ggplot2)
#bubble chart
ggplot(df_T23, 
       aes(x = category, y = ID, 
           size = p.adjust.log10,
           color = Regulation,
           shape=Regulation,
           #fill=Regulation
           fill=NES
       )
) +
  geom_point(
    stroke = 1.5,
    alpha = 0.9) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_manual(values = c("Up" = "black", "Down" = "black")) + # border color
  scale_shape_manual(values = c("Up" = 24, "Down" = 25)) + # shape legend
  scale_fill_gradient2(low = "#21318d", mid = "white", high = "red", midpoint = -1.3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.background = element_rect(fill = "white", color = "black"), 
    panel.border = element_rect(color = "black", fill=NA, size=1), 
    legend.key = element_blank(), 
    legend.key.size = unit(0.5, "cm"), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 14) 
  ) +
  labs(size = "-log10(FDR)", color = "Regulation", title = "Gene Set Enrichment Analysis")

ggsave("df_T23_g10.pdf", width = 9, height = 11,limitsize = FALSE)



# T2T1_g1
df_T2_g01 <- read.csv('GSEA.result.nodiff.new.20241224/5.diff_T2T1_g1_GSEA_symbols.csv',row.names = 1)
df_T3_g01 <- read.csv('GSEA.result.nodiff.new.20241224/7.diff_T3T1_g1_GSEA_symbols.csv',row.names = 1)

#筛选的通路
sx_T23_g1 <- c(
  'HALLMARK_MTORC1_SIGNALING',
  'HALLMARK_IL2_STAT5_SIGNALING',
  'REACTOME_INTERLEUKIN_3_INTERLEUKIN_5_AND_GM_CSF_SIGNALING',
  'REACTOME_INTERLEUKIN_7_SIGNALING',
  'REACTOME_ESR_MEDIATED_SIGNALING',
  'REACTOME_TCR_SIGNALING',
  'REACTOME_DECTIN_1_MEDIATED_NONCANONICAL_NF_KB_SIGNALING',
  'REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING',
  'REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT',
  'REACTOME_FC_EPSILON_RECEPTOR_FCERI_SIGNALING',
  'REACTOME_NEGATIVE_REGULATION_OF_NOTCH4_SIGNALING',
  'REACTOME_PCP_CE_PATHWAY',
  'REACTOME_INTERLEUKIN_RECEPTOR_SHC_SIGNALING',
  'REACTOME_SIGNALING_BY_NUCLEAR_RECEPTORS',
  'REACTOME_SIGNALING_BY_INTERLEUKINS',
  'REACTOME_DOWNSTREAM_SIGNALING_EVENTS_OF_B_CELL_RECEPTOR_BCR',
  'REACTOME_SIGNALING_BY_NOTCH',
  'HALLMARK_P53_PATHWAY',
  'REACTOME_CLEC7A_DECTIN_1_SIGNALING',
  'REACTOME_INTERLEUKIN_2_FAMILY_SIGNALING'
  
)

df_T2_g01 <- df_T2_g01[sx_T23_g1,]
df_T3_g01 <- df_T3_g01[sx_T23_g1,]
df_T2_g01 <- na.omit(df_T2_g01)
df_T3_g01 <- na.omit(df_T3_g01)

#T2
df_T2_g01_in <- df_T2_g01[,c('ID','NES','p.adjust')]
df_T2_g01_in$Regulation <- ifelse(df_T2_g01_in$NES>0,'Up','Down')
df_T2_g01_in$category <- 'T2'


#T3
df_T3_g01_in <- df_T3_g01[,c('ID','NES','p.adjust')]
df_T3_g01_in$Regulation <- ifelse(df_T3_g01_in$NES>0,'Up','Down')
df_T3_g01_in$category <- 'T3'




df_T23 <- rbind(df_T2_g01_in,df_T3_g01_in)
df_T23$p.adjust.log10 <- -log10(df_T23$p.adjust)
df_T23$ID <- tolower(df_T23$ID)

library(ggplot2)
# bubble plot
ggplot(df_T23, 
       aes(x = category, y = ID, 
           size = p.adjust.log10,
           color = Regulation,
           shape=Regulation,
           #fill=Regulation
           fill=NES
       )
) +
  geom_point(
    stroke = 1.5,
    alpha = 0.9) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_manual(values = c("Up" = "black", "Down" = "black")) + 
  scale_shape_manual(values = c("Up" = 24, "Down" = 25)) + 
  scale_fill_gradient2(low = "#21318d", mid = "white", high = "red", midpoint = -1.3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.background = element_rect(fill = "white", color = "black"), 
    panel.border = element_rect(color = "black", fill=NA, size=1), 
    legend.key = element_blank(), 
    legend.key.size = unit(0.5, "cm"), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 14)
  ) +
  labs(size = "-log10(FDR)", color = "Regulation", title = "Gene Set Enrichment Analysis")

ggsave("df_T23_g1.pdf", width = 9, height = 11,limitsize = FALSE)



