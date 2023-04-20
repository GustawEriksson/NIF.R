library("tximport")
library("readr")
library("tximportData")
library("pasilla")
library("DESeq2")
library("vsn")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("cowplot")
library("apeglm")
library("EnhancedVolcano")

# Setting parameters for the analysis
p.value.cutoff = 0.05
p.adj.cutoff = 0.05
log2FC.cutoff = 0.5
keep.samples = c("P14753_101", "P14753_102", "P14753_103", "P14753_104", "P14753_105", "P14753_106", # 3 weeks
                 "P14753_107", "P14753_108", "P14753_109", "P14753_110", "P14753_111", "P14753_112", # 6 weeks
                 "P14753_114", "P14753_115", "P14753_116", "P14753_118") # 18 weeks

# Generating the output directory
if (dir.exists(path = paste0("Output/1_DESeq2")) == FALSE) {
  print(paste0("Generating output directory Output/1_DESeq2"))
  dir.create(path = paste0("Output/1_DESeq2"), recursive = TRUE)
  Output.dir = paste0("Output/1_DESeq2/")
} else if (dir.exists(path = paste0("Output/1_DESeq2")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/1_DESeq2/")
} else {
  print("Error with output directory")
}

if (dir.exists(path = paste0("Output/2_Supplementary")) == FALSE) {
  print(paste0("Generating output directory Output/2_Supplementary"))
  dir.create(path = paste0("Output/2_Supplementary"), recursive = TRUE)
  Supp.dir = paste0("Output/2_Supplementary/")
} else if (dir.exists(path = paste0("Output/2_Supplementary")) == TRUE) {
  print("Directory exists")
  Supp.dir = paste0("Output/2_Supplementary/")
} else {
  print("Error with output directory")
}


# Read the data
cts = read.table(file = "Data/countsTabMatrix.txt", sep = "\t", header = TRUE)
rownames(cts) = cts[,1]
cts = cts[,-1]
coldata = openxlsx::read.xlsx(xlsxFile = "Data/Liver_experimental_design.xlsx", rowNames = TRUE)
coldata$Group <- factor(coldata$Group)

# Keep only liver samples
cts = cts[, c(keep.samples)]

# Generate DESeq2 dds object from matrix
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Group)
dds

# Run DESeq on the object
dds <- DESeq(dds)
res_w3 <- results(dds, contrast=c("Group", "NIF_W3", "Control_W3"))
res_w6 <- results(dds, contrast=c("Group", "NIF_W6", "Control_W6"))
res_w18 <- results(dds, contrast=c("Group", "NIF_W18", "Control_W18"))

# CHheck summaries of results
summary(res_w3)
summary(res_w6)
summary(res_w18)

# Do data transformation for plotting
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

meanSdPlot(assay(vsd)) # Better
meanSdPlot(assay(rld)) # Worse

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE)

# Quality control plotting
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Group)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup="Group")

# Initial plotting
plotMA(res_w3)

# Plotting the PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="Group")
