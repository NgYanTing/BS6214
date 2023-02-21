setwd("~/Documents/BS6214/Assignment1/")

library(DESeq2)
library(tximport)
library(ggplot2)
library(ComplexHeatmap)
library(ggrepel)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggvenn)
library(edgeR)

#### Prepare for analysis ####
# Make a phenotype table with 3 columns: sample    sex    condition
samplenames <- list.files("fastq_files")
samplenames <- gsub(".fastq.gz", "", samplenames)
pheno <- data.frame(sample = samplenames,
                    condition = c(rep("P19", 6), rep("CR", 6),
                                  rep("P5", 6), rep("P10", 6),
                                  rep("P15", 6)),
                    sex = rep(c("Male","Female"), 15)
                    )
write.csv(pheno, "phenotype.csv", row.names = F)

# Convert to factors and order accordingly
pheno$sample <- factor(pheno$sample)
pheno$condition <- factor(pheno$condition, levels = c("P19", "P15", "P10", "P5", "CR"))
pheno$sex <- factor(pheno$sex, levels = c("Male", "Female"))

# Get annotations
genes_to_transcripts <- read.delim("mus_musculus/transcripts_to_genes.txt",
                                   header = F)
colnames(genes_to_transcripts) <- c("TXNAME", "GENEID", "gene_name")

gtf <- read.delim("mus_musculus/Mus_musculus.GRCm38.96.gtf", header = F, skip = 5)[,c(3,9)]
gtf <- gtf[gtf$V3 == "gene",]
gtf$gene_name <- sapply(strsplit(gtf[[2]], ";", fixed = T), function(x) x[3])
gtf$gene_type <- sapply(strsplit(gtf[[2]], ";", fixed = T), function(x) x[5])
gtf$gene_name <- gsub("gene_name", "", gtf$gene_name)
gtf$gene_type <- gsub("gene_biotype", "", gtf$gene_type)
gtf$V3 <- NULL
gtf$V9 <- NULL
# Some columns have leading whitespace, so trim them off
gtf$gene_name <- trimws(gtf$gene_name)
gtf$gene_type <- trimws(gtf$gene_type)

# Append gene biotype to genes_to_transcripts
genes_to_transcripts <- merge(genes_to_transcripts, gtf, by = "gene_name")
rm(gtf)

# tximport
files <- file.path("./kallisto_out", pheno$sample, "abundance.h5")
names(files) <- paste0(pheno$sample)
# tx2gene should be a two-column data.frame linking transcript id
# (column 1) to gene id (column 2)
# This is genes_to_transcripts[,c(2,3)] in this case
txi.kallisto <- tximport(files, type = "kallisto",
                         txOut = F,
                         tx2gene = genes_to_transcripts[,c(2,3)])

#### PCA Plots ####
# Set up dds accounting for sex
dds <- DESeqDataSetFromTximport(txi.kallisto, pheno, ~sex + condition)
# Only keep protein_coding genes
pcg <- genes_to_transcripts[genes_to_transcripts$gene_type == "protein_coding",]$GENEID
dds <- dds[rownames(dds) %in% pcg,]
# Remove all lowly expressed genes (< 10 reads total)
dds <- dds[rowSums(counts(dds)) >= 10,]
# Remove all genes that are expressed only in less than 6 replicates
dds <- dds[rowSums(counts(dds) > 0) > 6, ]

## Plot PCA using vst
colorblind <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
vst <- vst(dds, blind = F)

png("pca.png", height = 400, width = 400, units = "px")
pcaData <- plotPCA(vst, intgroup = c("condition","sex"), returnData = T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$name <- gsub("_1", "", pcaData$name)
pcaData$label <- pcaData$name
pcaData$label <- c(pcaData$label[1], rep("",29))
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = sex)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = label)) +
  theme_bw(base_size = 18) +
  labs(title = "PCA Plot", color = "Condition", shape = "Sex",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = colorblind) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust= 0.5)
  ) + coord_fixed()
dev.off()
# There seems to be one outlier from the P19 group

# Remove outlier and re-setup dds
pheno <- pheno[-1,]
files <- files[-1]
txi.kallisto <- tximport(files, type = "kallisto",
                         txOut = F,
                         tx2gene = genes_to_transcripts[,c(2,3)])
dds <- DESeqDataSetFromTximport(txi.kallisto, pheno, ~sex + condition)
# Only keep protein_coding genes
pcg <- genes_to_transcripts[genes_to_transcripts$gene_type == "protein_coding",]$GENEID
dds <- dds[rownames(dds) %in% pcg,]
# Remove all lowly expressed genes (< 10 reads total)
dds <- dds[rowSums(counts(dds)) >= 10,]
# Remove all genes that are expressed only in less than 5 replicates
dds <- dds[rowSums(counts(dds) > 0) > 5, ]

# Re-plot PCA
vst <- vst(dds, blind = F)

png("pca2.png", height = 350, width = 400, units = "px")
pcaData <- plotPCA(vst, intgroup = c("condition","sex"), returnData = T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = sex)) +
  geom_point(size = 4) +
  theme_bw(base_size = 20) +
  labs(title = "PCA Plot", color = "Condition", shape = "Sex",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = colorblind) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust= 0.5)
  ) + coord_fixed()
dev.off()

#### DE gene analysis ####
## Calculate differential expression
X <- DESeq(dds)

# Results
process_results <- function(res) {
  resdf <- data.frame(res)
  resdf$GENEID <- rownames(resdf)
  row.names(resdf) <- 1:nrow(resdf)
  resdf <- merge(resdf, genes_to_transcripts, by = "GENEID")
  resdf <- resdf[,c("GENEID","gene_name","log2FoldChange","padj")]
  resdf <- unique(resdf)
  resdf <- resdf[order(resdf$padj, -abs(resdf$log2FoldChange)),]
}

# CR vs P19
res <- results(X, contrast=c("condition","CR","P19"), name = "test")
res <- lfcShrink(X, contrast=c("condition","CR","P19"), res = res, type = "ashr")
cr <- process_results(res)

# P5 vs P19
res <- results(X, contrast=c("condition","P5","P19"), name = "test")
res <- lfcShrink(X, contrast=c("condition","P5","P19"), res = res, type = "ashr")
p5 <- process_results(res)

# P10 vs P19
res <- results(X, contrast=c("condition","P10","P19"), name = "test")
res <- lfcShrink(X, contrast=c("condition","P10","P19"), res = res, type = "ashr")
p10 <- process_results(res)

# P15 vs P19
res <- results(X, contrast=c("condition","P15","P19"), name = "test")
res <- lfcShrink(X, contrast=c("condition","P15","P19"), res = res, type = "ashr")
p15 <- process_results(res)

## Assemble a table of top10 genes like the paper's table S2
top10 <- function(results, cond){
  tab <- results[,c("gene_name", "log2FoldChange", "padj")]
  colnames(tab) <- c(cond, "log2Fold", "padj")
  tab <- tab[order(tab$padj),]
  row.names(tab) <- 1:nrow(tab)
  tab <- tab[1:10,]
}
tabcr <- top10(cr, "CR")
tabp5 <- top10(p5, "5% Protein")
tabp10 <- top10(p10, "10% Protein")
tabp15 <- top10(p15, "15% Protein")
combined_table <- cbind(tabcr, tabp5, tabp10, tabp15)
write.csv(combined_table, "top10DEGs.csv", row.names = F)

#### GenAge ####
# Pro- and anti-longevity gene lists were downloaded from GenAge database

# Combine all results
process_each_res <- function(results, cond){
  tab <- results[,c("gene_name", "log2FoldChange","padj")]
  tab <- tab[,c("gene_name", "log2FoldChange")]
  # Scale values
  tab$log2FoldChange <- scale(tab$log2FoldChange)
  colnames(tab) <- c("gene_name", cond)
  return(tab)
}
rescr <- process_each_res(cr, "CR")
resp5 <- process_each_res(p5, "5%P")
resp10 <- process_each_res(p10, "10%P")
resp15 <- process_each_res(p15, "15%P")
all_res <- merge(resp15, resp10, by = "gene_name", all.x = T)
all_res <- merge(all_res, resp5, by = "gene_name", all.x = T)
all_res <- merge(all_res, rescr, by = "gene_name", all.x = T)
all_res <- unique(all_res)

# Anti-longevity genes in mice
anti <- read.delim("anti_longevity.tsv")[,1:3]
res_anti <- all_res[all_res$gene_name %in% anti$Gene.Symbol,]
hmap_mat <- as.matrix(res_anti[,-1])
row.names(hmap_mat) <- res_anti$gene_name

png("hmap_anti.png", height = 600, width = 400, units = "px")
Heatmap(hmap_mat, name = "Color Key", cluster_columns = F,
        clustering_method_rows = "complete", clustering_distance_rows = "euclidean",
        column_title = "Anti-longevity", column_title_gp = gpar(fontsize = 18),
        column_names_gp = gpar(fontsize = 16), column_names_rot = 0,
        column_names_centered = T,
        row_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Color Key",
                                    title_gp = gpar(size = 14)),
        col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
dev.off()

# Pro-longevity genes in mice
pro <- read.delim("pro_longevity.tsv")[,1:3]
res_pro <- all_res[all_res$gene_name %in% pro$Gene.Symbol,]
hmap_mat <- as.matrix(res_pro[,-1])
row.names(hmap_mat) <- res_pro$gene_name

png("hmap_pro.png", height = 600, width = 400, units = "px")
Heatmap(hmap_mat, name = "Color Key", cluster_columns = F,
        clustering_method_rows = "complete", clustering_distance_rows = "euclidean",
        column_title = "Pro-longevity", column_title_gp = gpar(fontsize = 18),
        column_names_gp = gpar(fontsize = 16), column_names_rot = 0,
        column_names_centered = T,
        row_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "Color Key",
                                    title_gp = gpar(size = 14)),
        col = circlize::colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")))
dev.off()

#### GSEA ####
clustPro <- function(results){
  tab <- results[,c("GENEID", "log2FoldChange", "padj")]
  # Remove .xx from the end of the ensembl IDs
  tab$GENEID <- gsub("\\..*","",tab$GENEID)
  # Order by smallest padj first
  # SetRank uses genes ranked by increasing p-value instead of log2FC
  # Clusterprofiler requires decreasing log2FC
  tab <- tab[order(-tab$log2FoldChange),]
  geneList <- tab$log2FoldChange
  names(geneList) <- tab$GENEID
  # omit NA values if any
  geneList <- na.omit(geneList)
  
  # GSEA GOBP
  gse <- gseGO(geneList = geneList, 
               ont = "ALL", 
               keyType = "ENSEMBL", 
               #nPerm = 10000, 
               minGSSize = 10, 
               maxGSSize = 1000, 
               pvalueCutoff = 0.05, 
               verbose = F, 
               OrgDb = org.Mm.eg.db, 
               pAdjustMethod = "BH")
  gse <- gse[,c("ID","Description", "pvalue", "p.adjust")]
  gse <- gse[gse$p.adjust < 0.05,]
  return(gse)
}

gse_cr <- clustPro(cr)
gse_p5 <- clustPro(p5)
gse_p10 <- clustPro(p10)
gse_p15 <- clustPro(p15)

# Save the top 24 entries; ID, description, adjusted p
write.csv(gse_cr[1:24,c("ID","Description", "p.adjust")], "gse_cr.csv", row.names = F)
write.csv(gse_p5[1:24,c("ID","Description", "p.adjust")], "gse_p5.csv", row.names = F)
write.csv(gse_p10[1:24,c("ID","Description", "p.adjust")], "gse_p10.csv", row.names = F)
write.csv(gse_p15[1:24,c("ID","Description", "p.adjust")], "gse_p15.csv", row.names = F)

#### Venn diagrams ####
# Trying to replicate a figure in the paper where the number of 
# common overrepresented GO categories for each experimental
# condition vs control are presented in a venn diagram

# Common GOs as venn diagram
x <- list(
  CR = gse_cr$ID,
  `5%P` = gse_p5$ID,
  `10%P` = gse_p10$ID,
  `15%P` = gse_p15$ID
)

png("venn.png", height = 400, width = 400, units = "px")
ggvenn(x,
       fill_color = c("indianred","royalblue","khaki1","palegreen"),
       stroke_size = 0.5,
       set_name_size = 8,
       text_size = 6,
       show_percentage = F
       )
dev.off()

#### Plot counts heatmap ####
# Should actually be using FPKM
pheno$cond_sex <- paste0(pheno$condition,pheno$sex)
pheno$cond_sex <- gsub("emale","",pheno$cond_sex)
pheno$cond_sex <- gsub("ale","",pheno$cond_sex)
pheno$cond_sex <- paste0(pheno$cond_sex, c(1,1,2,2,3,rep(c(1,1,2,2,3,3),4)))

cts <- counts(dds)
colnames(cts) <- pheno$cond_sex
cts <- t(scale(t(cts)))
cts <- cts[,c("P19F1","P19F2","P19F3","P19M1","P19M2",
              "P15F1","P15F2","P15F3","P15M1","P15M2","P15M3",
              "P10F1","P10F2","P10F3","P10M1","P10M2","P10M3",
              "P5F1","P5F2","P5F3","P5M1","P5M2","P5M3",
              "CRF1","CRF2","CRF3","CRM1","CRM2","CRM3")]
# Filter counts to keep only those that are significantly DE from P19 control
de_genes <- cr[cr$padj < 0.05,]$GENEID
de_genes <- c(de_genes, p5[p5$padj < 0.05,]$GENEID)
de_genes <- c(de_genes, p10[p10$padj < 0.05,]$GENEID)
de_genes <- c(de_genes, p15[p15$padj < 0.05,]$GENEID)
de_genes <- unique(de_genes)
cts <- cts[rownames(cts) %in% de_genes,]

png("hmap_all.png", height = 550, width = 550, units = "px")
Heatmap(cts, name = "Color Key", cluster_columns = F,
        show_row_names = F,
        heatmap_legend_param = list(title = "Color Key",
                                    title_gp = gpar(size = 14)),
        clustering_method_rows = "complete", clustering_distance_rows = "euclidean",
        col = circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
        )
dev.off()

#### edgeR ####
# Use edgeR for comparison to the DEGs identified from DESeq2
cts <- txi.kallisto$counts
normMat <- txi.kallisto$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)

# filtering using the design information
design <- model.matrix(~sex + condition, data = pheno)
keep <- filterByExpr(y, design)
y <- y[keep, ]

# Drop all non-protein-coding genes
pcg <- genes_to_transcripts[genes_to_transcripts$gene_type == "protein_coding",]$GENEID
y <- y[rownames(y) %in% pcg,]

y <- estimateDisp(y)
group <- pheno$condition
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
P19vsP15 <- glmQLFTest(fit, coef=2)
P19vsP10 <- glmQLFTest(fit, coef=3)
P19vsP5 <- glmQLFTest(fit, coef=4)
P19vsCR <- glmQLFTest(fit, coef=5)

top10 <- function(df, cond){
  df$GENEID <- row.names(df)
  df <- df[,c("GENEID", "logFC", "FDR")]
  df <- merge(df, genes_to_transcripts[,c(1,3)], by = "GENEID")
  df <- unique(df[,c("gene_name", "logFC", "FDR")])
  colnames(df) <- c(cond, "log2Fold", "padj")
  df <- df[order(df$padj, -df$log2Fold),]
  row.names(df) <- 1:nrow(df)
  return(df)
}

tabcr <- top10(as.data.frame(topTags(P19vsCR)), "CR")
tabp5 <- top10(as.data.frame(topTags(P19vsP5)), "5% Protein")
tabp10 <- top10(as.data.frame(topTags(P19vsP10)), "10% Protein")
tabp15 <- top10(as.data.frame(topTags(P19vsP15)), "15% Protein")
combined_table <- cbind(tabcr, tabp5, tabp10, tabp15)
write.csv(combined_table, "top10DEGs_edgeR.csv", row.names = F)

sessionInfo()
