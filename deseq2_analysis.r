library(DESeq2)
counts_data <- read.csv('C:/Users/user/Downloads/2117ensemble.csv', row.names = 1)
colData <- read.csv("C:/users/user/OneDrive/Documents/2117sample_infor.csv")
row.names(colData) = colData$Tags
all(colnames(counts_data) %in% rownames(colData))
all(colnames(counts_data) == rownames(colData))
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~ disease)
results(dds, contrast = c("disease", "NGT", "T2D"))
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$disease <- relevel(dds$disease, ref = "T2D")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
plotMA(dds)re
write.csv(df, "DESeqResults.csv", row.names = FALSE)
transformed_counts <- counts(dds,normalized = TRUE)
write.csv( transformed_counts, "transformed_counts.csv")


#volcano plot
library(EnhancedVolcano)
EnhancedVolcano(res,
        lab = rownames(res),
        x = 'log2FoldChange',
        y = 'padj',
        xlab = 'log2 fold change',
        ylab = '-log10 adjusted p-value',
        title = 'T2D versus NGT',
        legendLabels = c('NS', 'log2 FC', 'adj p-value', 'adj p-value and log2 FC'),
        pCutoffCol = 'padj',
        pCutoff=0.05,
        FCcutoff = 1,
        
        pointSize = 3.0,
        labSize = 6.0)

#GO enrichment
 go_bp <- enrichGO(
     gene = degs$Tags,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",  
     pvalueCutoff = 0.05,   
      ont = "BP"    

 go_cc <- enrichGO(
     gene = degs$Tags,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",   
     pvalueCutoff = 0.05,    
     qvalueCutoff = 0.2,
      ont = "CC"    
 )

 go_mf <- enrichGO(
     gene = degs$Tags,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",   
     pvalueCutoff = 0.05,    
     qvalueCutoff = 0.2,
      ont = "MF"  
 )
dotplot(go_cc, showCategory=20)






