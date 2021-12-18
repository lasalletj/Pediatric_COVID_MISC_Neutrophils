Pediatric\_COVID\_MISC\_Neutrophils
================
Tom LaSalle

This document contains all the code necessary to generate the plots for
all figures pertaining to neutrophil RNA-seq data (Figure 2, Figure 3,
Supplementary Figure S1, Supplementary Figure S2, Supplementary Figure
S3). Plots are subsequently edited in Adobe Illustrator to produce the
final figures.

Load the necessary libraries:

Set the color palette:

``` r
well <- "#008300"
covid <- "#007fc7"
misc <- "#d00022"
vermillion <- rgb(213,94,0,max=255) 
bluishgreen <- rgb(0,158,115,max=255)
yellow <- rgb(240,228,66,max=255)
blue <- rgb(0,114,178,max=255)
orange <- rgb(230,159,0,max=255)
skyblue <- rgb(86,180,233,max=255)
lightgray <- rgb(211,211,211,max=255)
pink <- "#c64dd1"
```

Import the metadata and gene expression matrices. We remove duplicated
samples and samples that failed quality control from the metadata table:

``` r
prefix <- "~/Downloads/Pedi_COVID_MISC_Supplementary_Tables/"
qc <- read.xlsx(paste0(prefix,"data-file-S1.xlsx"), sheet = 2)
qc <- qc[complete.cases(qc$Public.Sample.ID),]
cibersortx <- read.xlsx(paste0(prefix,"data-file-S1.xlsx"), sheet = 3)
cibersortx <- cibersortx[complete.cases(cibersortx$Public.Sample.ID),]
metadata <- cbind(cibersortx,qc[,2:ncol(qc)])
remove(qc,cibersortx)
metadata <- metadata[-grep("A",metadata$Public.Sample.ID),] # Remove duplicated samples
```

Import RNA-seq Count and TPM matrices, filtering out lowly expressed
genes:

``` r
Count <- read.xlsx(paste0(prefix,"data-file-S3.xlsx"), sheet = 2)
Count <- Count[-grep("PAR_Y",Count$gene_id),]
rownames(Count) <- gsub("\\..*","",Count$gene_id)
Count <- Count[,-c(1,2)]

TPM <- read.xlsx(paste0(prefix,"data-file-S3.xlsx"), sheet = 3)
TPM <- TPM[-grep("PAR_Y",TPM$gene_id),]
rownames(TPM) <- gsub("\\..*","",TPM$gene_id)
TPM <- TPM[,-c(1,2)]

genepc <- read.delim(paste0(prefix,"Ensembl_to_Symbol.txt"))
logTPM <- log2(TPM + 1)

tf <- rowSums(TPM > 0.1) > ncol(TPM)*.2
TPM <- TPM[tf,]
Count <- Count[tf,]
logTPM <- logTPM[tf,]
tf <- rowSums(Count >= 6) > ncol(Count)*.2
TPM <- TPM[tf,]
Count <- Count[tf,]
logTPM <- logTPM[tf,]
```

First we check that our bulk neutrophil RNA-seq samples are high purity
using a digital cytometry method titled CIBERSORTx. The count matrix is
uploaded to the online tool at <https://cibersortx.stanford.edu/> along
with the cell type signature matrix. We read in the results here:

``` r
# Transform this data in %
cibersort_temp <- metadata[,colnames(metadata) %in% c("disease","Mature_Neutrophil","Immature_Neutrophil","Monocyte","T_NK","B","Plasmablast","Total_Neutrophil")]
cibersort_temp$disease <- factor(cibersort_temp$disease, levels = c("HC","COVID","MISC"))
cibersort_temp <- cibersort_temp[rownames(cibersort_temp) %in% rownames(metadata),]
cibersort_sorted <- cibersort_temp[order(cibersort_temp$disease,cibersort_temp$Total_Neutrophil),]
cibersort_sorted <- cibersort_sorted[,-c(1,4)]
data_percentage <- t(cibersort_sorted*100)
coul <- c("forestgreen","tomato","skyblue","slateblue3","gray","black")
```

**Figure S1A:**

``` r
barplot(data_percentage, col=coul , border=NA, xlab="group")
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

We represent the same data as boxplots:

``` r
p1 <- ggplot(metadata, aes(x = factor(disease, levels = c("HC","COVID","MISC")), y = Total_Neutrophil*100, fill = factor(disease, levels = c("HC","COVID","MISC")))) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, alpha = 0.3, pch = 16) + stat_compare_means(label.y = 2) + theme_bw() + xlab("") + ylab("CIBERSORTx Percentage (%)") + ggtitle("Total Neutrophil") + theme(panel.grid = element_blank(), legend.position = "none") + scale_fill_manual(values = c(well,covid,misc)) + scale_y_continuous(limits = c(0,100.1))
p1$labels$fill <- ""

p2 <- ggplot(metadata, aes(x = factor(disease, levels = c("HC","COVID","MISC")), y = Mature_Neutrophil*100, fill = factor(disease, levels = c("HC","COVID","MISC")))) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, alpha = 0.3, pch = 16) + stat_compare_means(label.y = 2) + theme_bw() + xlab("") + ylab("") + ggtitle("Mature Neutrophil") + theme(panel.grid = element_blank(), legend.position = "none") + scale_fill_manual(values = c(well,covid,misc)) + scale_y_continuous(limits = c(0,100.1))
p2$labels$fill <- ""

p3 <- ggplot(metadata, aes(x = factor(disease, levels = c("HC","COVID","MISC")), y = Immature_Neutrophil*100, fill = factor(disease, levels = c("HC","COVID","MISC")))) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, alpha = 0.3, pch = 16) + stat_compare_means(label.y = 2) + theme_bw() + xlab("") + ylab("") + ggtitle("Immature Neutrophil") + theme(panel.grid = element_blank(), legend.position = "none") + scale_fill_manual(values = c(well,covid,misc)) + scale_y_continuous(limits = c(0,100.1))
p3$labels$fill <- ""

p4 <- ggplot(metadata, aes(x = factor(disease, levels = c("HC","COVID","MISC")), y = Monocyte*100, fill = factor(disease, levels = c("HC","COVID","MISC")))) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, alpha = 0.3, pch = 16) + stat_compare_means(label.y = 2) + theme_bw() + xlab("") + ylab("CIBERSORTx Percentage (%)") + ggtitle("Monocyte") + theme(panel.grid = element_blank(), legend.position = "none") + scale_fill_manual(values = c(well,covid,misc)) + scale_y_continuous(limits = c(0,100.1))
p4$labels$fill <- ""

p5 <- ggplot(metadata, aes(x = factor(disease, levels = c("HC","COVID","MISC")), y = T_NK*100, fill = factor(disease, levels = c("HC","COVID","MISC")))) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, alpha = 0.3, pch = 16) + stat_compare_means(label.y = 2) + theme_bw() + xlab("") + ylab("") + ggtitle("T_NK") + theme(panel.grid = element_blank(), legend.position = "none") + scale_fill_manual(values = c(well,covid,misc)) + scale_y_continuous(limits = c(0,100.1))
p5$labels$fill <- ""

p6 <- ggplot(metadata, aes(x = factor(disease, levels = c("HC","COVID","MISC")), y = B*100, fill = factor(disease, levels = c("HC","COVID","MISC")))) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, alpha = 0.3, pch = 16) + stat_compare_means(label.y = 2) + theme_bw() + xlab("") + ylab("") + ggtitle("B") + theme(panel.grid = element_blank(), legend.position = "none") + scale_fill_manual(values = c(well,covid,misc)) + scale_y_continuous(limits = c(0,100.1))
p6$labels$fill <- ""

p7 <- ggplot(metadata, aes(x = factor(disease, levels = c("HC","COVID","MISC")), y = Plasmablast*100, fill = factor(disease, levels = c("HC","COVID","MISC")))) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, alpha = 0.3, pch = 16) + stat_compare_means(label.y = 2) + theme_bw() + xlab("") + ylab("") + ggtitle("Plasmablast") + theme(panel.grid = element_blank(), legend.position = "none") + scale_fill_manual(values = c(well,covid,misc)) + scale_y_continuous(limits = c(0,100.1))
p7$labels$fill <- ""
```

**Figure S1B:**

``` r
plot_grid(p1,p2,p3,NA,p4,p5,p6,p7,ncol=4)
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

We also confirmed high purity of samples using flow cytometry in
**Figure S1C**. Next we generate a UMAP of all RNA-seq samples that
passed quality control.

``` r
#set.seed(10101)
#pcg.umap <- umap(t(logTPM))
#metadata$umap1 <- pcg.umap$layout[,1]
#metadata$umap2 <- pcg.umap$layout[,2]
umap_coordinates <- read.xlsx(paste0(prefix,"umap_coordinates.xlsx"))
metadata <- cbind(metadata,umap_coordinates[,2:3])

p1 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = factor(disease, levels = c("HC","COVID","MISC")))) + geom_point(size = 2.5) + theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank()) + ggtitle("UMAP") + xlab("") + ylab("") + scale_colour_manual(values = c(well,covid,misc)) + coord_fixed(ratio = 2)
p1$labels$colour <- "Diagnosis"
```

**Figure S1D:**

``` r
p1
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

We can superimpose the CIBERSORTx fractions to see the groupings of
mature and immature neutrophils.

``` r
myPalette <- colorRampPalette((brewer.pal(9, "RdYlBu")))
sc <- scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(0,100))
p1 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = Mature_Neutrophil*100)) + sc + geom_point(size = 2.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), legend.position = "none") + ggtitle("Mature Neutrophil") + coord_fixed(ratio = .95) + xlab("") + ylab("") + coord_fixed(ratio = 2)
p1$labels$colour <- ""

p2 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = Immature_Neutrophil*100)) + sc + geom_point(size = 2.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), legend.position = "none") + ggtitle("Immature Neutrophil") + coord_fixed(ratio = .95) + xlab("") + ylab("") + coord_fixed(ratio = 2)
p2$labels$colour <- ""

p3 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = T_NK*100)) + sc + geom_point(size = 2.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), legend.position = "none") + ggtitle("T/NK") + coord_fixed(ratio = .95) + xlab("") + ylab("") + coord_fixed(ratio = 2)
p3$labels$colour <- ""

p4 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = Monocyte*100)) + sc + geom_point(size = 2.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), legend.position = "none") + ggtitle("Monocyte") + coord_fixed(ratio = .95) + xlab("") + ylab("") + coord_fixed(ratio = 2)
p4$labels$colour <- ""

p5 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = B*100)) + sc + geom_point(size = 2.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), legend.position = "none") + ggtitle("B") + coord_fixed(ratio = .95) + xlab("") + ylab("") + coord_fixed(ratio = 2)
p5$labels$colour <- ""

p6 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = Plasmablast*100)) + sc + geom_point(size = 2.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), legend.position = "none") + ggtitle("Plasmablast") + coord_fixed(ratio = .95) + xlab("") + ylab("") + coord_fixed(ratio = 2)
p6$labels$colour <- ""
```

**Figure S1E:**

``` r
plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Now we begin our analysis of neutrophils in acute pediatric COVID-19 by
performing a differential expression analysis and GSEA against healthy
controls.

``` r
coldata <- metadata
coldata <- coldata[coldata$disease %in% c("HC","COVID"),]
coldata$disease <- factor(coldata$disease, levels = c("HC","COVID"))
Count_select <- Count[,which(colnames(Count) %in% coldata$Public.Sample.ID)]
logTPM_select <- logTPM[,which(colnames(logTPM) %in% coldata$Public.Sample.ID)]
Count_rounded <- round(Count_select)

dds <- DESeqDataSetFromMatrix(countData = Count_rounded, colData = coldata, design = ~ disease)

dds <- DESeq(dds)
res <- results(dds, name="disease_COVID_vs_HC")
filenam <- "COVID_vs_HC"
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")
temp <- genepc[which(genepc$Gene.stable.ID %in% rownames(res)),]
res$gene_symbol <- matrix(0L, nrow = nrow(res))
for (i in 1:nrow(res)){
  if (rownames(res)[i] %in% temp$Gene.stable.ID){
    res$gene_symbol[i] <- temp$Gene.name[which(rownames(res)[i] == temp$Gene.stable.ID)]
  } else {
    res$gene_symbol[i] <- rownames(res)[i]
  }
}
res$rank <- sign(res$log2FoldChange)*(-1)*log10(res$pvalue)
res <- res[complete.cases(res),]
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- res[rev(order(res$rank)),]

resgsea <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")

ranking <- resgsea[,"rank"]
names(ranking) <- resgsea$gene_symbol

allgenesetgmt.file <- c(gmtPathways(paste0(prefix,"all_gene_sets.gmt")))
set.seed(15001)
allgenesetfgseaRes <- fgsea(pathways = allgenesetgmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
allgenesetfgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(allgenesetfgseaRes), ncol = 1)
for (i in 1:nrow(allgenesetfgseaRes)){
  allgenesetfgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(allgenesetfgseaRes$leadingEdge[i]), collapse = ", "))
}
allgenesetfgseaRes <- allgenesetfgseaRes[rev(order(allgenesetfgseaRes$NES)),]

neustategmt.file <- c(gmtPathways(paste0(prefix,"neutrophil_state_gene_sets.gmt")))
set.seed(15001)
neustatefgseaRes <- fgsea(pathways = neustategmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
neustatefgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(neustatefgseaRes), ncol = 1)
for (i in 1:nrow(neustatefgseaRes)){
  neustatefgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(neustatefgseaRes$leadingEdge[i]), collapse = ", "))
}
neustatefgseaRes <- neustatefgseaRes[rev(order(neustatefgseaRes$NES)),]

sheets <- list("DESeq2" = cbind(rownames(res),res[,c(7,1,2,3,4,5,6,8)]), "PATHWAYS_GSEA" = allgenesetfgseaRes[,c(1,2,3,4,5,6,7,9)], "NEU_STATE_GSEA" = neustatefgseaRes[,c(1,2,3,4,5,6,7,9)]) 
write_xlsx(sheets, paste0(prefix,filenam,".xlsx"))

resordered <- res[order(res$rank),]

log2fc <- as.numeric(resordered$log2FoldChange)
log10p <- as.numeric(-1*log10(resordered$pvalue))
padj <- resordered$padj
rank <- resordered$rank
symbol <- resordered$gene_symbol
combo <- cbind(log2fc,log10p,padj,symbol,rank)
rownames(combo) <- rownames(resordered)
colnames(combo) <- c("log2fc","log10p","padj","gene_symbol","rank")
combo <- as.data.frame(combo)
combo$rank <- as.numeric(combo$rank)
combo$log10p <- as.numeric(combo$log10p)
combo$log2fc <- as.numeric(combo$log2fc)
combo$significance <- as.numeric(combo$padj < 0.05)
combo$significance <- as.factor(combo$significance)

combo$labels <- 0
combo$labels[combo$gene_symbol %in% c("KLHDC3","GLIPR1","REX1BD","NDUFA13","SCAP","CCL2","CCRL2","MOB3C","SLC43A3","ELMO2","GPR84","CCL4L2","TNFSF8","SERPINB9","MX1","IRF7","IFIT3","TCN2","AP3B2")] <- 1
combo$labels <- as.factor(combo$labels)

options(ggrepel.max.overlaps = Inf)
my.cols <- c(well, lightgray, covid)
plot1 <- ggplot(combo, aes(x = log2fc, y = log10p)) + geom_point(data = subset(combo, significance == 1 & log2fc < 0), colour = my.cols[1]) + geom_point(data = subset(combo, significance == 0), colour = my.cols[2]) + geom_point(data = subset(combo, significance == 1 & log2fc > 0), colour = my.cols[3]) + geom_text_repel(data = subset(combo, labels == 1), aes(label = as.character(gene_symbol)))  + theme_bw() + ylab("-Log10(p-value)") + xlab("Log2(FC)") + ggtitle("COVID vs. Healthy") + annotate("text", x=1.3, y=0, label= "COVID", colour = my.cols[3]) + annotate("text", x=-1.2, y=0, label= "Healthy", colour = my.cols[1]) + coord_fixed(ratio = .6, xlim = c(-1.5,3.5)) + theme(panel.grid = element_blank())
```

**Figure 2A:**

``` r
plot1
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
allgenesetfgseaResordered <- allgenesetfgseaRes[(order(allgenesetfgseaRes$NES)),]
allgenesetfgseaResordered <- allgenesetfgseaResordered[allgenesetfgseaResordered$pathway %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB","ARDS_UP_JUSS","GO_DEFENSE_RESPONSE_TO_VIRUS","GO_OXIDATIVE_PHOSPHORYLATION","ARDS_DOWN_JUSS","GO_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN","GO_TRANSLATIONAL_INITIATION","GO_INNER_MITOCHONDRIAL_MEMBRANE_ORGANIZATION"),]

allgenesetfgseaResordered$idx <- 1:nrow(allgenesetfgseaResordered)
allgenesetfgseaResordered$sign <- allgenesetfgseaResordered$NES > 0
allgenesetfgseaResordered$sign <- mapvalues(allgenesetfgseaResordered$sign, from = c("FALSE","TRUE"), to = c("Healthy","COVID"))
```

**Figure 2B:**

``` r
ggplot(allgenesetfgseaResordered, aes(x = factor(idx), y = NES)) + geom_col(aes(fill = sign)) + coord_flip() + theme_bw() + theme(panel.grid = element_blank()) + scale_fill_manual(values = c(covid, well))
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Next we want to compare the pediatric neutrophil profiles to
previously-defined neutrophil states in adults in the context of
SARS-CoV-2 infection. We start by introducing the cohort from LaSalle et
al. 2021 (biorXiv) utilizing publicly available data from Supplementary
Table 1 presented in a new plot.

``` r
adult_neu <- read.xlsx(paste0(prefix,"LaSalle_Clusters.xlsx"))
adult_neu$cluster_neuhi <- mapvalues(adult_neu$cluster_neuhi, from = c("1","2","3","4","5","6","7"), to = c("Pro-Neu","NFkB+","PD-L1+ ISG+","Immature","G-MDSC-like","ISG+","Neu-Lo"))
p1 <- ggplot(data = adult_neu, aes(x = umap1, y = umap2, colour = factor(cluster_neuhi, levels = c("Pro-Neu","NFkB+","PD-L1+ ISG+","Immature","G-MDSC-like","ISG+","Neu-Lo")))) + geom_point(size = 1) + theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank()) + ggtitle("UMAP") + xlab("") + ylab("") + scale_colour_manual(values = c(orange,skyblue,bluishgreen,yellow,blue,vermillion,pink)) + coord_fixed(ratio = .75)
p1$labels$colour <- ""
```

**Figure 2C:**

``` r
p1
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Next we can assign each sample a score according to the marker genes of
the various subtypes (excluding Neu-Lo, which was not a cluster but a
designation of lower purity).

``` r
subtypescore <- function(subtype){
  subtypegenes <- neustategmt.file[[subtype]]
  ids <- matrix(0L, nrow = length(subtypegenes), ncol = 1)
  subtypegenes <- as.data.frame(cbind(subtypegenes, ids))
  colnames(subtypegenes) <- c("symbol","id")
  temp <- genepc[genepc$Gene.stable.ID %in% rownames(logTPM),]
  for (i in 1:nrow(subtypegenes)){
    symbol <- subtypegenes$symbol[i]
    ids <- temp[temp$Gene.name == symbol,]
    if (nrow(ids)==1){
      subtypegenes$id[i] <- ids$Gene.stable.ID[1]
    }
  }
  subtypegenes <- subtypegenes[subtypegenes$id %in% rownames(logTPM),]

  logTPM_select <- logTPM
  set.seed(150001)
  geneset <- subtypegenes$id
  SC <- matrix(0L, nrow = 1, ncol = ncol(logTPM_select))
  colnames(SC) <- colnames(logTPM_select)

  normdata <- logTPM_select
  normdata.t <- t(logTPM_select)
  normdata <- as.data.frame(normdata)
  normdata.t <- as.data.frame(normdata.t)

  Ea <- as.data.frame(rowSums(logTPM_select))
  colnames(Ea) <- "Ea"
  Ea$bin <- matrix(0L, ncol = 1, nrow = nrow(Ea))
  Ea <- Ea[order(Ea$Ea),]
  Ea$bin <- c(rep(1,617),rep(2,616),rep(3,616),rep(4,616),rep(5,617),rep(6,616),rep(7,616),rep(8,616),rep(9,616),rep(10,617),rep(11,616),rep(12,616),rep(13,616),rep(14,616),rep(15,617),rep(16,616),rep(17,616),rep(18,616),rep(19,616),rep(20,617),rep(21,616),rep(22,616),rep(23,616),rep(24,616),rep(25,617))

  geneset <- geneset[which(geneset %in% colnames(normdata.t))]
  selectedgenes.t <- subset(normdata.t, select = geneset)
  selectedgenes <- as.data.frame(t(selectedgenes.t))

  binnum <- Ea[geneset[1],2]
  selectedbin <- Ea[Ea$bin == binnum,]
  selection <- sample.int(nrow(selectedbin),100,replace = FALSE)
  controlgenelist <- rownames(selectedbin[selection,])
  binnumbers <- binnum
  
  for (j in 2:length(geneset)){
    binnum <- Ea[geneset[j],2]
    selectedbin <- Ea[Ea$bin == binnum,]
    selection <- sample.int(nrow(selectedbin),100, replace = FALSE)
    controlgenes <- rownames(selectedbin[selection,])
    controlgenelist <- c(controlgenelist,controlgenes)
    binnumbers <- paste(binnumbers,binnum)
  }
  
  controlgenelist <- unique(controlgenelist)
  controlgenes.t <- subset(normdata.t, select = controlgenelist)
  controlgenes <- as.data.frame(t(controlgenes.t))
 
  i = 1
  j = 1
  for (i in 1:length(SC)){
  
    ErGji <- sum(selectedgenes[,i])/length(selectedgenes[,i])
  
    ErGjconti <- sum(controlgenes[,i])/length(controlgenes[,i])
  
    SC[i] <- ErGji - ErGjconti
  }
  return(t(SC))
}

metadata$NMF1 <- subtypescore("NMF1")
metadata$NMF2 <- subtypescore("NMF2")
metadata$NMF3 <- subtypescore("NMF3")
metadata$NMF4 <- subtypescore("NMF4")
metadata$NMF5 <- subtypescore("NMF5")
metadata$NMF6 <- subtypescore("NMF6")

myPalette <- colorRampPalette((brewer.pal(9, "RdBu")))
sc <- scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(min(metadata$NMF1,metadata$NMF2,metadata$NMF3,metadata$NMF4,metadata$NMF5,metadata$NMF6),max(metadata$NMF1,metadata$NMF2,metadata$NMF3,metadata$NMF4,metadata$NMF5,metadata$NMF6)))

metadata$NMF1 <- as.numeric(metadata$NMF1)
p1 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = as.numeric(NMF1))) + sc + geom_point(size = 3.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), axis.line = element_blank(), legend.position = "none") + ggtitle("NMF1") + coord_fixed(ratio = 1.9) + 
  xlab("") + ylab("")
p1$labels$colour <- ""

metadata$NMF2 <- as.numeric(metadata$NMF2)
p2 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = as.numeric(NMF2))) + sc + geom_point(size = 3.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), axis.line = element_blank(), legend.position = "none") + ggtitle("NMF2") + coord_fixed(ratio = 1.9) + 
  xlab("") + ylab("") 
p2$labels$colour <- ""

metadata$NMF3 <- as.numeric(metadata$NMF3)
p3 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = as.numeric(NMF3))) + sc + geom_point(size = 3.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), axis.line = element_blank(), legend.position = "none") + ggtitle("NMF3") + coord_fixed(ratio = 1.9) + 
  xlab("") + ylab("")
p3$labels$colour <- ""

metadata$NMF4 <- as.numeric(metadata$NMF4)
p4 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = as.numeric(NMF4))) + sc + geom_point(size = 3.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), axis.line = element_blank(), legend.position = "none") + ggtitle("NMF4") + coord_fixed(ratio = 1.9) + 
  xlab("") + ylab("") 
p4$labels$colour <- ""

metadata$NMF5 <- as.numeric(metadata$NMF5)
p5 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = as.numeric(NMF5))) + sc + geom_point(size = 3.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), axis.line = element_blank(), legend.position = "none") + ggtitle("NMF5") + coord_fixed(ratio = 1.9) + 
  xlab("") + ylab("") 
p5$labels$colour <- ""

metadata$NMF6 <- as.numeric(metadata$NMF6)
p6 <- ggplot(data = metadata, aes(x = umap1, y = umap2, colour = as.numeric(NMF6))) + sc + geom_point(size = 3.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank(), panel.grid = element_blank(), axis.line = element_blank(), legend.position = "none") + ggtitle("NMF6") + coord_fixed(ratio = 1.9) + 
  xlab("") + ylab("")  
p6$labels$colour <- ""
```

**Figure 2D:**

``` r
plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

We then look for enrichment of those neutrophil states in our
differential expression result between acute pediatric COVID-19 and
healthy controls.

``` r
getEnrichmentDataframe <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2) {
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    combo <- as.data.frame(cbind(tops,bottoms))
    combo$average <- matrix(0L, nrow = nrow(combo), ncol = 1)
    for (p in 1:nrow(combo)){
      combo$average[p] <- (combo$tops[p]+combo$bottoms[p])/2
    }
    combo <- cbind(combo,pathway)
    return(combo)
}

res <- res[rev(order(res$rank)),]
ranking <- res[,"rank"]
names(ranking) <- res$gene_symbol
NMFscorepathways <- c("NMF1","NMF3","NMF5","NMF6")

df <- getEnrichmentDataframe(neustategmt.file[[NMFscorepathways[1]]], ranking)
NMFscoredf <- cbind(df$pathway,df$average,rep(NMFscorepathways[1],nrow(df)))
colnames(NMFscoredf) <- c("rank","enrichment","pathway")

for (i in 2:length(NMFscorepathways)){
  df <- getEnrichmentDataframe(neustategmt.file[[NMFscorepathways[i]]],
               ranking)
  temp <- cbind(df$pathway,df$average,rep(NMFscorepathways[i],nrow(df)))
  colnames(temp) <- c("rank","enrichment","pathway")
  NMFscoredf <- as.data.frame(rbind(NMFscoredf,temp))
}

my.cols <- c(orange,bluishgreen,blue,vermillion)
p1 <- ggplot(as.data.frame(NMFscoredf), aes(x = as.numeric(rank), y = as.numeric(enrichment), colour = pathway)) + geom_point(aes(shape = pathway)) + theme_bw() + scale_colour_manual(values = my.cols) + theme(panel.grid.minor = element_blank(), legend.text=element_text(size=6)) + xlab("Rank in Gene List") + ylab("Running Enrichment Score") + ggtitle("NMF Scores") + scale_shape_manual(values = c(16,15,17,18))
```

**Figure 2E:**

``` r
p1
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

Next we move on to MIS-C versus healthy control comparisons. We start
with differential gene expression and GSEA.

``` r
coldata <- metadata
coldata <- coldata[coldata$disease %in% c("HC","MISC"),]
coldata$disease <- factor(coldata$disease, levels = c("HC","MISC"))
Count_select <- Count[,which(colnames(Count) %in% coldata$Public.Sample.ID)]
logTPM_select <- logTPM[,which(colnames(logTPM) %in% coldata$Public.Sample.ID)]
Count_rounded <- round(Count_select)

dds <- DESeqDataSetFromMatrix(countData = Count_rounded, colData = coldata, design = ~ disease)

dds <- DESeq(dds)
res <- results(dds, name="disease_MISC_vs_HC")
filenam <- "MISC_vs_HC"
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")
temp <- genepc[which(genepc$Gene.stable.ID %in% rownames(res)),]
res$gene_symbol <- matrix(0L, nrow = nrow(res))
for (i in 1:nrow(res)){
  if (rownames(res)[i] %in% temp$Gene.stable.ID){
    res$gene_symbol[i] <- temp$Gene.name[which(rownames(res)[i] == temp$Gene.stable.ID)]
  } else {
    res$gene_symbol[i] <- rownames(res)[i]
  }
}
res$rank <- sign(res$log2FoldChange)*(-1)*log10(res$pvalue)
res <- res[complete.cases(res),]
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- res[rev(order(res$rank)),]

resgsea <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")

ranking <- resgsea[,"rank"]
names(ranking) <- resgsea$gene_symbol

allgenesetgmt.file <- c(gmtPathways(paste0(prefix,"all_gene_sets.gmt")))
set.seed(15001)
allgenesetfgseaRes <- fgsea(pathways = allgenesetgmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
allgenesetfgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(allgenesetfgseaRes), ncol = 1)
for (i in 1:nrow(allgenesetfgseaRes)){
  allgenesetfgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(allgenesetfgseaRes$leadingEdge[i]), collapse = ", "))
}
allgenesetfgseaRes <- allgenesetfgseaRes[rev(order(allgenesetfgseaRes$NES)),]

neustategmt.file <- c(gmtPathways(paste0(prefix,"neutrophil_state_gene_sets.gmt")))
set.seed(15001)
neustatefgseaRes <- fgsea(pathways = neustategmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
neustatefgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(neustatefgseaRes), ncol = 1)
for (i in 1:nrow(neustatefgseaRes)){
  neustatefgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(neustatefgseaRes$leadingEdge[i]), collapse = ", "))
}
neustatefgseaRes <- neustatefgseaRes[rev(order(neustatefgseaRes$NES)),]

sheets <- list("DESeq2" = cbind(rownames(res),res[,c(7,1,2,3,4,5,6,8)]), "PATHWAYS_GSEA" = allgenesetfgseaRes[,c(1,2,3,4,5,6,7,9)], "NEU_STATE_GSEA" = neustatefgseaRes[,c(1,2,3,4,5,6,7,9)]) 
write_xlsx(sheets, paste0(prefix,filenam,".xlsx"))

resordered <- res[order(res$rank),]

log2fc <- as.numeric(resordered$log2FoldChange)
log10p <- as.numeric(-1*log10(resordered$pvalue))
padj <- resordered$padj
rank <- resordered$rank
symbol <- resordered$gene_symbol
combo <- cbind(log2fc,log10p,padj,symbol,rank)
rownames(combo) <- rownames(resordered)
colnames(combo) <- c("log2fc","log10p","padj","gene_symbol","rank")
combo <- as.data.frame(combo)
combo$rank <- as.numeric(combo$rank)
combo$log10p <- as.numeric(combo$log10p)
combo$log2fc <- as.numeric(combo$log2fc)
combo$significance <- as.numeric(combo$padj < 0.05)
combo$significance <- as.factor(combo$significance)

combo$labels <- 0
combo$labels[combo$gene_symbol %in% c("PLEKHO1","SNRNP70","DTX4","DEFA1","GPR84","ACER3","LDHA","EXOSC4","FLVCR2","ENO1","C3AR1","HILPDA","CDKN1A","FKBP4")] <- 1
combo$labels <- as.factor(combo$labels)

options(ggrepel.max.overlaps = Inf)
my.cols <- c(well, lightgray, misc)
plot1 <- ggplot(combo, aes(x = log2fc, y = log10p)) + geom_point(data = subset(combo, significance == 1 & log2fc < 0), colour = my.cols[1]) + geom_point(data = subset(combo, significance == 0), colour = my.cols[2]) + geom_point(data = subset(combo, significance == 1 & log2fc > 0), colour = my.cols[3]) + geom_text_repel(data = subset(combo, labels == 1), aes(label = as.character(gene_symbol)))  + theme_bw() + ylab("-Log10(p-value)") + xlab("Log2(FC)") + ggtitle("MIS-C vs. Healthy") + annotate("text", x=1.3, y=0, label= "MIS-C", colour = my.cols[3]) + annotate("text", x=-1.2, y=0, label= "Healthy", colour = my.cols[1]) + coord_fixed(ratio = .9, xlim = c(-1.2,5)) + theme(panel.grid = element_blank())
```

**Figure 3A:**

``` r
plot1
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
allgenesetfgseaResordered <- allgenesetfgseaRes[(order(allgenesetfgseaRes$NES)),]
allgenesetfgseaResordered <- allgenesetfgseaResordered[allgenesetfgseaResordered$pathway %in% c("GO_CHROMOSOME_ORGANIZATION","ARDS_DOWN_JUSS","GO_POSITIVE_REGULATION_OF_GTPASE_ACTIVITY","GO_CHROMATIN_ORGANIZATION","ARDS_UP_JUSS","HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_MTORC1_SIGNALING","REACTOME_NEUTROPHIL_DEGRANULATION","HALLMARK_GLYCOLYSIS","GO_RESPONSE_TO_INTERLEUKIN_1","HALLMARK_TNFA_SIGNALING_VIA_NFKB"),]

allgenesetfgseaResordered$idx <- 1:nrow(allgenesetfgseaResordered)
allgenesetfgseaResordered$sign <- allgenesetfgseaResordered$NES > 0
allgenesetfgseaResordered$sign <- mapvalues(allgenesetfgseaResordered$sign, from = c("FALSE","TRUE"), to = c("Healthy","COVID"))
```

**Figure 3B:**

``` r
ggplot(allgenesetfgseaResordered, aes(x = factor(idx), y = NES)) + geom_col(aes(fill = sign)) + coord_flip() + theme_bw() + theme(panel.grid = element_blank()) + scale_fill_manual(values = c(misc, well))
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

We then look for enrichment of neutrophil states in our differential
expression result between MIS-C and healthy controls.

``` r
res <- res[rev(order(res$rank)),]
ranking <- res[,"rank"]
names(ranking) <- res$gene_symbol
NMFscorepathways <- c("NMF4","NMF5")

df <- getEnrichmentDataframe(neustategmt.file[[NMFscorepathways[1]]], ranking)
NMFscoredf <- cbind(df$pathway,df$average,rep(NMFscorepathways[1],nrow(df)))
colnames(NMFscoredf) <- c("rank","enrichment","pathway")

for (i in 2:length(NMFscorepathways)){
  df <- getEnrichmentDataframe(neustategmt.file[[NMFscorepathways[i]]],
               ranking)
  temp <- cbind(df$pathway,df$average,rep(NMFscorepathways[i],nrow(df)))
  colnames(temp) <- c("rank","enrichment","pathway")
  NMFscoredf <- as.data.frame(rbind(NMFscoredf,temp))
}

my.cols <- c(yellow,blue)
p1 <- ggplot(as.data.frame(NMFscoredf), aes(x = as.numeric(rank), y = as.numeric(enrichment), colour = pathway)) + geom_point(aes(shape = pathway)) + theme_bw() + scale_colour_manual(values = my.cols) + theme(panel.grid.minor = element_blank(), legend.text=element_text(size=6)) + xlab("Rank in Gene List") + ylab("Running Enrichment Score") + ggtitle("NMF Scores") + scale_shape_manual(values = c(16,15))
```

**Figure 3C:**

``` r
p1
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

To confirm our findings, we perform a direct comparison of MIS-C and
acute pediatric COVID-19 neutrophils.

``` r
coldata <- metadata
coldata <- coldata[coldata$disease %in% c("COVID","MISC"),]
coldata$disease <- factor(coldata$disease, levels = c("COVID","MISC"))
Count_select <- Count[,which(colnames(Count) %in% coldata$Public.Sample.ID)]
logTPM_select <- logTPM[,which(colnames(logTPM) %in% coldata$Public.Sample.ID)]
Count_rounded <- round(Count_select)

dds <- DESeqDataSetFromMatrix(countData = Count_rounded, colData = coldata, design = ~ disease)

dds <- DESeq(dds)
res <- results(dds, name="disease_MISC_vs_COVID")
filenam <- "MISC_vs_COVID"
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")
temp <- genepc[which(genepc$Gene.stable.ID %in% rownames(res)),]
res$gene_symbol <- matrix(0L, nrow = nrow(res))
for (i in 1:nrow(res)){
  if (rownames(res)[i] %in% temp$Gene.stable.ID){
    res$gene_symbol[i] <- temp$Gene.name[which(rownames(res)[i] == temp$Gene.stable.ID)]
  } else {
    res$gene_symbol[i] <- rownames(res)[i]
  }
}
res$rank <- sign(res$log2FoldChange)*(-1)*log10(res$pvalue)
res <- res[complete.cases(res),]
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- res[rev(order(res$rank)),]

resgsea <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")

ranking <- resgsea[,"rank"]
names(ranking) <- resgsea$gene_symbol

allgenesetgmt.file <- c(gmtPathways(paste0(prefix,"all_gene_sets.gmt")))
set.seed(15001)
allgenesetfgseaRes <- fgsea(pathways = allgenesetgmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
allgenesetfgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(allgenesetfgseaRes), ncol = 1)
for (i in 1:nrow(allgenesetfgseaRes)){
  allgenesetfgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(allgenesetfgseaRes$leadingEdge[i]), collapse = ", "))
}
allgenesetfgseaRes <- allgenesetfgseaRes[rev(order(allgenesetfgseaRes$NES)),]

neustategmt.file <- c(gmtPathways(paste0(prefix,"neutrophil_state_gene_sets.gmt")))
set.seed(15001)
neustatefgseaRes <- fgsea(pathways = neustategmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
neustatefgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(neustatefgseaRes), ncol = 1)
for (i in 1:nrow(neustatefgseaRes)){
  neustatefgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(neustatefgseaRes$leadingEdge[i]), collapse = ", "))
}
neustatefgseaRes <- neustatefgseaRes[rev(order(neustatefgseaRes$NES)),]

sheets <- list("DESeq2" = cbind(rownames(res),res[,c(7,1,2,3,4,5,6,8)]), "PATHWAYS_GSEA" = allgenesetfgseaRes[,c(1,2,3,4,5,6,7,9)], "NEU_STATE_GSEA" = neustatefgseaRes[,c(1,2,3,4,5,6,7,9)]) 
write_xlsx(sheets, paste0(prefix,filenam,".xlsx"))

resordered <- res[order(res$rank),]

log2fc <- as.numeric(resordered$log2FoldChange)
log10p <- as.numeric(-1*log10(resordered$pvalue))
padj <- resordered$padj
rank <- resordered$rank
symbol <- resordered$gene_symbol
combo <- cbind(log2fc,log10p,padj,symbol,rank)
rownames(combo) <- rownames(resordered)
colnames(combo) <- c("log2fc","log10p","padj","gene_symbol","rank")
combo <- as.data.frame(combo)
combo$rank <- as.numeric(combo$rank)
combo$log10p <- as.numeric(combo$log10p)
combo$log2fc <- as.numeric(combo$log2fc)
combo$significance <- as.numeric(combo$padj < 0.05)
combo$significance <- as.factor(combo$significance)

combo$labels <- 0
combo$labels[combo$gene_symbol %in% c("CCL2","OTOF","MX1","FCGR2C","OASL","IFIT2","IFIT5","FCGR2B","MX2","IFI44L","CCL4L2","DDX58","TNFSF10","IFI44","IDH3G","AP1M1","AP2S1","TUFM","HDAC3","PDHA1","SLC25A5")] <- 1
combo$labels <- as.factor(combo$labels)

options(ggrepel.max.overlaps = Inf)
my.cols <- c(covid, lightgray, misc)
plot1 <- ggplot(combo, aes(x = log2fc, y = log10p)) + geom_point(data = subset(combo, significance == 1 & log2fc < 0), colour = my.cols[1]) + geom_point(data = subset(combo, significance == 0), colour = my.cols[2]) + geom_point(data = subset(combo, significance == 1 & log2fc > 0), colour = my.cols[3]) + geom_text_repel(data = subset(combo, labels == 1), aes(label = as.character(gene_symbol)))  + theme_bw() + ylab("-Log10(p-value)") + xlab("Log2(FC)") + ggtitle("MIS-C vs. COVID-19") + annotate("text", x=1.3, y=0, label= "MIS-C", colour = my.cols[3]) + annotate("text", x=-1.2, y=0, label= "COVID-19", colour = my.cols[1]) + coord_fixed(ratio = .9, xlim = c(-6.5,6.5)) + theme(panel.grid = element_blank())
```

**Figure S3A:**

``` r
plot1
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

We confirm enrichment of the previously-identified neutrophil states:

``` r
res <- res[rev(order(res$rank)),]
ranking <- res[,"rank"]
names(ranking) <- res$gene_symbol
NMFscorepathways <- c("NMF1","NMF3","NMF4","NMF5","NMF6")

df <- getEnrichmentDataframe(neustategmt.file[[NMFscorepathways[1]]], ranking)
NMFscoredf <- cbind(df$pathway,df$average,rep(NMFscorepathways[1],nrow(df)))
colnames(NMFscoredf) <- c("rank","enrichment","pathway")

for (i in 2:length(NMFscorepathways)){
  df <- getEnrichmentDataframe(neustategmt.file[[NMFscorepathways[i]]],
               ranking)
  temp <- cbind(df$pathway,df$average,rep(NMFscorepathways[i],nrow(df)))
  colnames(temp) <- c("rank","enrichment","pathway")
  NMFscoredf <- as.data.frame(rbind(NMFscoredf,temp))
}

my.cols <- c(orange,bluishgreen,yellow,blue,vermillion)
p1 <- ggplot(as.data.frame(NMFscoredf), aes(x = as.numeric(rank), y = as.numeric(enrichment), colour = pathway)) + geom_point(aes(shape = pathway)) + theme_bw() + scale_colour_manual(values = my.cols) + theme(panel.grid.minor = element_blank(), legend.text=element_text(size=6)) + xlab("Rank in Gene List") + ylab("Running Enrichment Score") + ggtitle("NMF Scores") + scale_shape_manual(values = c(16,15,17,18,16))
```

**Figure S3B:**

``` r
p1
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

We are also interested in identifying unique marker genes that
distinguish acute pediatric COVID-19, MIS-C, and healthy controls. To do
so we perform three differential expression analyses: (MIS-C and
COVID-19) vs. Healthy, (MIS-C and Healthy) vs. COVID-19, and (COVID-19
and Healthy) vs. MIS-C.

``` r
metadata$COVID_group <- as.numeric(metadata$disease == "COVID")
metadata$COVID_group <- mapvalues(metadata$COVID_group, from = c("0","1"), to = c("other","COVID"))
metadata$COVID_group <- factor(metadata$COVID_group, levels = c("other","COVID"))
metadata$MISC_group <- as.numeric(metadata$disease == "MISC")
metadata$MISC_group <- mapvalues(metadata$MISC_group, from = c("0","1"), to = c("other","MISC"))
metadata$MISC_group <- factor(metadata$MISC_group, levels = c("other","MISC"))
metadata$HC_group <- as.numeric(metadata$disease == "HC")
metadata$HC_group <- mapvalues(metadata$HC_group, from = c("0","1"), to = c("other","HC"))
metadata$HC_group <- factor(metadata$HC_group, levels = c("other","HC"))
```

``` r
coldata <- metadata
Count_select <- Count[,which(colnames(Count) %in% coldata$Public.Sample.ID)]
logTPM_select <- logTPM[,which(colnames(logTPM) %in% coldata$Public.Sample.ID)]
Count_rounded <- round(Count_select)

dds <- DESeqDataSetFromMatrix(countData = Count_rounded, colData = coldata, design = ~ COVID_group)

dds <- DESeq(dds)
res <- results(dds, name="COVID_group_COVID_vs_other")
filenam <- "COVID_vs_other"
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")
temp <- genepc[which(genepc$Gene.stable.ID %in% rownames(res)),]
res$gene_symbol <- matrix(0L, nrow = nrow(res))
for (i in 1:nrow(res)){
  if (rownames(res)[i] %in% temp$Gene.stable.ID){
    res$gene_symbol[i] <- temp$Gene.name[which(rownames(res)[i] == temp$Gene.stable.ID)]
  } else {
    res$gene_symbol[i] <- rownames(res)[i]
  }
}
res$rank <- sign(res$log2FoldChange)*(-1)*log10(res$pvalue)
res <- res[complete.cases(res),]
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- res[rev(order(res$rank)),]

resgsea <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")

ranking <- resgsea[,"rank"]
names(ranking) <- resgsea$gene_symbol

allgenesetgmt.file <- c(gmtPathways(paste0(prefix,"all_gene_sets.gmt")))
set.seed(15001)
allgenesetfgseaRes <- fgsea(pathways = allgenesetgmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
allgenesetfgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(allgenesetfgseaRes), ncol = 1)
for (i in 1:nrow(allgenesetfgseaRes)){
  allgenesetfgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(allgenesetfgseaRes$leadingEdge[i]), collapse = ", "))
}
allgenesetfgseaRes <- allgenesetfgseaRes[rev(order(allgenesetfgseaRes$NES)),]

neustategmt.file <- c(gmtPathways(paste0(prefix,"neutrophil_state_gene_sets.gmt")))
set.seed(15001)
neustatefgseaRes <- fgsea(pathways = neustategmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
neustatefgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(neustatefgseaRes), ncol = 1)
for (i in 1:nrow(neustatefgseaRes)){
  neustatefgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(neustatefgseaRes$leadingEdge[i]), collapse = ", "))
}
neustatefgseaRes <- neustatefgseaRes[rev(order(neustatefgseaRes$NES)),]

sheets <- list("DESeq2" = cbind(rownames(res),res[,c(7,1,2,3,4,5,6,8)]), "PATHWAYS_GSEA" = allgenesetfgseaRes[,c(1,2,3,4,5,6,7,9)], "NEU_STATE_GSEA" = neustatefgseaRes[,c(1,2,3,4,5,6,7,9)]) 
write_xlsx(sheets, paste0(prefix,filenam,".xlsx"))

coldata <- metadata
Count_select <- Count[,which(colnames(Count) %in% coldata$Public.Sample.ID)]
logTPM_select <- logTPM[,which(colnames(logTPM) %in% coldata$Public.Sample.ID)]
Count_rounded <- round(Count_select)

dds <- DESeqDataSetFromMatrix(countData = Count_rounded, colData = coldata, design = ~ MISC_group)

dds <- DESeq(dds)
res <- results(dds, name="MISC_group_MISC_vs_other")
filenam <- "MISC_vs_other"
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")
temp <- genepc[which(genepc$Gene.stable.ID %in% rownames(res)),]
res$gene_symbol <- matrix(0L, nrow = nrow(res))
for (i in 1:nrow(res)){
  if (rownames(res)[i] %in% temp$Gene.stable.ID){
    res$gene_symbol[i] <- temp$Gene.name[which(rownames(res)[i] == temp$Gene.stable.ID)]
  } else {
    res$gene_symbol[i] <- rownames(res)[i]
  }
}
res$rank <- sign(res$log2FoldChange)*(-1)*log10(res$pvalue)
res <- res[complete.cases(res),]
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- res[rev(order(res$rank)),]

resgsea <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")

ranking <- resgsea[,"rank"]
names(ranking) <- resgsea$gene_symbol

allgenesetgmt.file <- c(gmtPathways(paste0(prefix,"all_gene_sets.gmt")))
set.seed(15001)
allgenesetfgseaRes <- fgsea(pathways = allgenesetgmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
allgenesetfgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(allgenesetfgseaRes), ncol = 1)
for (i in 1:nrow(allgenesetfgseaRes)){
  allgenesetfgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(allgenesetfgseaRes$leadingEdge[i]), collapse = ", "))
}
allgenesetfgseaRes <- allgenesetfgseaRes[rev(order(allgenesetfgseaRes$NES)),]

neustategmt.file <- c(gmtPathways(paste0(prefix,"neutrophil_state_gene_sets.gmt")))
set.seed(15001)
neustatefgseaRes <- fgsea(pathways = neustategmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
neustatefgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(neustatefgseaRes), ncol = 1)
for (i in 1:nrow(neustatefgseaRes)){
  neustatefgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(neustatefgseaRes$leadingEdge[i]), collapse = ", "))
}
neustatefgseaRes <- neustatefgseaRes[rev(order(neustatefgseaRes$NES)),]

sheets <- list("DESeq2" = cbind(rownames(res),res[,c(7,1,2,3,4,5,6,8)]), "PATHWAYS_GSEA" = allgenesetfgseaRes[,c(1,2,3,4,5,6,7,9)], "NEU_STATE_GSEA" = neustatefgseaRes[,c(1,2,3,4,5,6,7,9)]) 
write_xlsx(sheets, paste0(prefix,filenam,".xlsx"))

coldata <- metadata
Count_select <- Count[,which(colnames(Count) %in% coldata$Public.Sample.ID)]
logTPM_select <- logTPM[,which(colnames(logTPM) %in% coldata$Public.Sample.ID)]
Count_rounded <- round(Count_select)

dds <- DESeqDataSetFromMatrix(countData = Count_rounded, colData = coldata, design = ~ HC_group)

dds <- DESeq(dds)
res <- results(dds, name="HC_group_HC_vs_other")
filenam <- "HC_vs_other"
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")
temp <- genepc[which(genepc$Gene.stable.ID %in% rownames(res)),]
res$gene_symbol <- matrix(0L, nrow = nrow(res))
for (i in 1:nrow(res)){
  if (rownames(res)[i] %in% temp$Gene.stable.ID){
    res$gene_symbol[i] <- temp$Gene.name[which(rownames(res)[i] == temp$Gene.stable.ID)]
  } else {
    res$gene_symbol[i] <- rownames(res)[i]
  }
}
res$rank <- sign(res$log2FoldChange)*(-1)*log10(res$pvalue)
res <- res[complete.cases(res),]
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- res[rev(order(res$rank)),]

resgsea <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")

ranking <- resgsea[,"rank"]
names(ranking) <- resgsea$gene_symbol

allgenesetgmt.file <- c(gmtPathways(paste0(prefix,"all_gene_sets.gmt")))
set.seed(15001)
allgenesetfgseaRes <- fgsea(pathways = allgenesetgmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
allgenesetfgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(allgenesetfgseaRes), ncol = 1)
for (i in 1:nrow(allgenesetfgseaRes)){
  allgenesetfgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(allgenesetfgseaRes$leadingEdge[i]), collapse = ", "))
}
allgenesetfgseaRes <- allgenesetfgseaRes[rev(order(allgenesetfgseaRes$NES)),]

neustategmt.file <- c(gmtPathways(paste0(prefix,"neutrophil_state_gene_sets.gmt")))
set.seed(15001)
neustatefgseaRes <- fgsea(pathways = neustategmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
neustatefgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(neustatefgseaRes), ncol = 1)
for (i in 1:nrow(neustatefgseaRes)){
  neustatefgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(neustatefgseaRes$leadingEdge[i]), collapse = ", "))
}
neustatefgseaRes <- neustatefgseaRes[rev(order(neustatefgseaRes$NES)),]

sheets <- list("DESeq2" = cbind(rownames(res),res[,c(7,1,2,3,4,5,6,8)]), "PATHWAYS_GSEA" = allgenesetfgseaRes[,c(1,2,3,4,5,6,7,9)], "NEU_STATE_GSEA" = neustatefgseaRes[,c(1,2,3,4,5,6,7,9)]) 
write_xlsx(sheets, paste0(prefix,filenam,".xlsx"))
```

We now build the heatmap of the marker genes. We add clinical metadata
to the top of each column. The time from treatment until sample
collection is not public information; all other clinical metadata is
found in Supplementary Tables 1 and 2. In this heatmap, “Treatment”
refers to supplemental O2 for acute pediatric COVID-19 patients and IVIG
for MIS-C patients. Further rearrangement of column order was performed
in Adobe Illustrator.

**Figure S3:**

``` r
pheatmap(logTPM_temp, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = TRUE, scale = "row", color=colfunc(length(breaksList)), breaks = breaksList, annotation_col = coldata[,colnames(coldata) %in% c("disease","Treatment","Cardiac")], fontsize_row = 3, cellwidth = 6, cellheight = 3, border_color = NA)
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

The final RNA-seq comparison is between MIS-C patients with cardiac
involvement and MIS-C patients without cardiac involvement.

``` r
coldata <- metadata
coldata <- coldata[coldata$disease %in% c("MISC"),]
coldata$Cardiac <- factor(coldata$Cardiac, levels = c("no","yes"))
Count_select <- Count[,which(colnames(Count) %in% coldata$Public.Sample.ID)]
logTPM_select <- logTPM[,which(colnames(logTPM) %in% coldata$Public.Sample.ID)]
Count_rounded <- round(Count_select)

dds <- DESeqDataSetFromMatrix(countData = Count_rounded, colData = coldata, design = ~ Cardiac)

dds <- DESeq(dds)
res <- results(dds, name="Cardiac_yes_vs_no")
filenam <- "Cardiac_yes_vs_no"
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")
temp <- genepc[which(genepc$Gene.stable.ID %in% rownames(res)),]
res$gene_symbol <- matrix(0L, nrow = nrow(res))
for (i in 1:nrow(res)){
  if (rownames(res)[i] %in% temp$Gene.stable.ID){
    res$gene_symbol[i] <- temp$Gene.name[which(rownames(res)[i] == temp$Gene.stable.ID)]
  } else {
    res$gene_symbol[i] <- rownames(res)[i]
  }
}
res$rank <- sign(res$log2FoldChange)*(-1)*log10(res$pvalue)
res <- res[complete.cases(res),]
write.table(res,paste0(prefix,filenam,".txt"),sep = "\t")
res <- res[rev(order(res$rank)),]

resgsea <- read.delim(paste0(prefix,filenam,".txt"),sep = "\t")

ranking <- resgsea[,"rank"]
names(ranking) <- resgsea$gene_symbol

allgenesetgmt.file <- c(gmtPathways(paste0(prefix,"all_gene_sets.gmt")))
set.seed(15001)
allgenesetfgseaRes <- fgsea(pathways = allgenesetgmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
allgenesetfgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(allgenesetfgseaRes), ncol = 1)
for (i in 1:nrow(allgenesetfgseaRes)){
  allgenesetfgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(allgenesetfgseaRes$leadingEdge[i]), collapse = ", "))
}
allgenesetfgseaRes <- allgenesetfgseaRes[rev(order(allgenesetfgseaRes$NES)),]

neustategmt.file <- c(gmtPathways(paste0(prefix,"neutrophil_state_gene_sets.gmt")))
set.seed(15001)
neustatefgseaRes <- fgsea(pathways = neustategmt.file, 
                  stats = ranking,
                  minSize=1,
                  maxSize=5000,
                  eps = 0)
neustatefgseaRes$LeadingEdgeGenes <- matrix(0L, nrow = nrow(neustatefgseaRes), ncol = 1)
for (i in 1:nrow(neustatefgseaRes)){
  neustatefgseaRes$LeadingEdgeGenes[i] <- as.character(paste(unlist(neustatefgseaRes$leadingEdge[i]), collapse = ", "))
}
neustatefgseaRes <- neustatefgseaRes[rev(order(neustatefgseaRes$NES)),]

sheets <- list("DESeq2" = cbind(rownames(res),res[,c(7,1,2,3,4,5,6,8)]), "PATHWAYS_GSEA" = allgenesetfgseaRes[,c(1,2,3,4,5,6,7,9)], "NEU_STATE_GSEA" = neustatefgseaRes[,c(1,2,3,4,5,6,7,9)]) 
write_xlsx(sheets, paste0(prefix,filenam,".xlsx"))

resordered <- res[order(res$rank),]


res <- readxl::read_xlsx(paste0(prefix,"Cardiac_Yes_vs_No.xlsx"), sheet = 1)
resordered <- res[rev(order(res$rank)),]
resordered <- resordered[resordered$padj < 0.05,]
resordered$idx <- rev(1:nrow(resordered))
resordered$sign <- resordered$log2FoldChange > 0
resordered$sign <- mapvalues(resordered$sign, from = c("FALSE","TRUE"), to = c("No","Yes"))
p1 <- ggplot(resordered, aes(x = factor(idx), y = rank)) + geom_col(aes(fill = sign)) + geom_text(aes(label = gene_symbol)) + coord_flip() + theme_bw() + theme(panel.grid = element_blank()) + scale_fill_manual(values = c(yellow,orange))
p1$labels$fill <- "Cardiac"
```

**Figure 3D:**

``` r
p1
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

Finally we highlight a few important genes.

``` r
metadata_temp <- metadata[metadata$disease %in% c("MISC"),]

genepc <- genepc[genepc$Gene.stable.ID %in% rownames(logTPM),]

goi <- "HLA-DQB1"
a <- genepc$Gene.stable.ID[which(genepc$Gene.name == goi)][1]
logTPM_temp <- logTPM[,colnames(logTPM) %in% metadata_temp$Public.Sample.ID]
metadata_temp$GOI <- t(logTPM_temp[which(rownames(logTPM_temp) == a),])
p1 <- ggplot(metadata_temp, aes(x = Cardiac, y = GOI, fill = Cardiac)) + theme_bw() + theme(legend.position = "none", panel.grid = element_blank()) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, alpha = 0.3, pch = 16) + xlab("") + ylab("") + ggtitle(goi) + scale_fill_manual(values = c(yellow,orange)) + coord_fixed(ratio = .32)

goi <- "HLA-DRA"
a <- genepc$Gene.stable.ID[which(genepc$Gene.name == goi)][1]
logTPM_temp <- logTPM[,colnames(logTPM) %in% metadata_temp$Public.Sample.ID]
metadata_temp$GOI <- t(logTPM_temp[which(rownames(logTPM_temp) == a),])
p2 <- ggplot(metadata_temp, aes(x = Cardiac, y = GOI, fill = Cardiac)) + theme_bw() + theme(legend.position = "none", panel.grid = element_blank()) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, alpha = 0.3, pch = 16) + xlab("") + ylab("") + ggtitle(goi) + scale_fill_manual(values = c(yellow,orange)) + coord_fixed(ratio = .32)

goi <- "HLA-DMB"
a <- genepc$Gene.stable.ID[which(genepc$Gene.name == goi)][1]
logTPM_temp <- logTPM[,colnames(logTPM) %in% metadata_temp$Public.Sample.ID]
metadata_temp$GOI <- t(logTPM_temp[which(rownames(logTPM_temp) == a),])
p3 <- ggplot(metadata_temp, aes(x = Cardiac, y = GOI, fill = Cardiac)) + theme_bw() + theme(legend.position = "none", panel.grid = element_blank()) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, alpha = 0.3, pch = 16) + xlab("") + ylab("") + ggtitle(goi) + scale_fill_manual(values = c(yellow,orange)) + coord_fixed(ratio = .38)

goi <- "CD163"
a <- genepc$Gene.stable.ID[which(genepc$Gene.name == goi)][1]
logTPM_temp <- logTPM[,colnames(logTPM) %in% metadata_temp$Public.Sample.ID]
metadata_temp$GOI <- t(logTPM_temp[which(rownames(logTPM_temp) == a),])
p4 <- ggplot(metadata_temp, aes(x = Cardiac, y = GOI, fill = Cardiac)) + theme_bw() + theme(legend.position = "none", panel.grid = element_blank()) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, alpha = 0.3, pch = 16) + xlab("") + ylab("") + ggtitle(goi) + scale_fill_manual(values = c(yellow,orange)) + coord_fixed(ratio = .28)

goi <- "ZBTB16"
a <- genepc$Gene.stable.ID[which(genepc$Gene.name == goi)][1]
logTPM_temp <- logTPM[,colnames(logTPM) %in% metadata_temp$Public.Sample.ID]
metadata_temp$GOI <- t(logTPM_temp[which(rownames(logTPM_temp) == a),])
p5 <- ggplot(metadata_temp, aes(x = Cardiac, y = GOI, fill = Cardiac)) + theme_bw() + theme(legend.position = "none", panel.grid = element_blank()) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, alpha = 0.3, pch = 16) + xlab("") + ylab("") + ggtitle(goi) + scale_fill_manual(values = c(yellow,orange)) + coord_fixed(ratio = .295)

goi <- "OLAH"
a <- genepc$Gene.stable.ID[which(genepc$Gene.name == goi)][1]
logTPM_temp <- logTPM[,colnames(logTPM) %in% metadata_temp$Public.Sample.ID]
metadata_temp$GOI <- t(logTPM_temp[which(rownames(logTPM_temp) == a),])
p6 <- ggplot(metadata_temp, aes(x = Cardiac, y = GOI, fill = Cardiac)) + theme_bw() + theme(legend.position = "none", panel.grid = element_blank()) + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, alpha = 0.3, pch = 16) + xlab("") + ylab("") + ggtitle(goi) + scale_fill_manual(values = c(yellow,orange)) + coord_fixed(ratio = .26)
```

**Figure 3E:**

``` r
plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)
```

![](Pediatric_COVID_MISC_Neutrophils_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] fgsea_1.18.0                DESeq2_1.32.0              
    ##  [3] SummarizedExperiment_1.22.0 Biobase_2.52.0             
    ##  [5] MatrixGenerics_1.4.3        matrixStats_0.60.1         
    ##  [7] GenomicRanges_1.44.0        GenomeInfoDb_1.28.2        
    ##  [9] IRanges_2.26.0              S4Vectors_0.30.0           
    ## [11] purrr_0.3.4                 writexl_1.4.0              
    ## [13] openxlsx_4.2.4              ggrepel_0.9.1              
    ## [15] umap_0.2.7.0                cowplot_1.1.1              
    ## [17] Pigengene_1.18.0            BiocStyle_2.20.2           
    ## [19] graph_1.70.0                BiocGenerics_0.38.0        
    ## [21] pheatmap_1.0.12             dplyr_1.0.7                
    ## [23] plyr_1.8.6                  RColorBrewer_1.1-2         
    ## [25] stringr_1.4.0               knitr_1.33                 
    ## [27] ggpubr_0.4.0                gridExtra_2.3              
    ## [29] ggplot2_3.3.5              
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1           backports_1.2.1        fastmatch_1.1-3       
    ##   [4] Hmisc_4.5-0            splines_4.1.1          BiocParallel_1.26.2   
    ##   [7] digest_0.6.27          foreach_1.5.1          htmltools_0.5.2       
    ##  [10] GO.db_3.13.0           gdata_2.18.0           fansi_0.5.0           
    ##  [13] magrittr_2.0.1         checkmate_2.0.0        memoise_2.0.0         
    ##  [16] cluster_2.1.2          doParallel_1.0.16      annotate_1.70.0       
    ##  [19] fastcluster_1.2.3      Biostrings_2.60.2      askpass_1.1           
    ##  [22] jpeg_0.1-9             colorspace_2.0-2       blob_1.2.2            
    ##  [25] haven_2.4.3            xfun_0.25              crayon_1.4.1          
    ##  [28] RCurl_1.98-1.4         jsonlite_1.7.2         libcoin_1.0-8         
    ##  [31] genefilter_1.74.0      impute_1.66.0          survival_3.2-13       
    ##  [34] iterators_1.0.13       glue_1.4.2             gtable_0.3.0          
    ##  [37] zlibbioc_1.38.0        XVector_0.32.0         DelayedArray_0.18.0   
    ##  [40] car_3.0-11             Rgraphviz_2.36.0       abind_1.4-5           
    ##  [43] scales_1.1.1           mvtnorm_1.1-2          DBI_1.1.1             
    ##  [46] rstatix_0.7.0          Rcpp_1.0.7             xtable_1.8-4          
    ##  [49] htmlTable_2.2.1        Cubist_0.3.0           reticulate_1.20       
    ##  [52] foreign_0.8-81         bit_4.0.4              preprocessCore_1.54.0 
    ##  [55] Formula_1.2-4          htmlwidgets_1.5.3      httr_1.4.2            
    ##  [58] ellipsis_0.3.2         farver_2.1.0           XML_3.99-0.7          
    ##  [61] pkgconfig_2.0.3        nnet_7.3-16            locfit_1.5-9.4        
    ##  [64] utf8_1.2.2             dynamicTreeCut_1.63-1  labeling_0.4.2        
    ##  [67] tidyselect_1.1.1       rlang_0.4.11           reshape2_1.4.4        
    ##  [70] AnnotationDbi_1.54.1   munsell_0.5.0          cellranger_1.1.0      
    ##  [73] tools_4.1.1            cachem_1.0.6           generics_0.1.0        
    ##  [76] RSQLite_2.2.8          broom_0.7.9            evaluate_0.14         
    ##  [79] fastmap_1.1.0          yaml_2.2.1             bit64_4.0.5           
    ##  [82] zip_2.2.0              KEGGREST_1.32.0        compiler_4.1.1        
    ##  [85] rstudioapi_0.13        bnlearn_4.6.1          curl_4.3.2            
    ##  [88] png_0.1-7              ggsignif_0.6.2         geneplotter_1.70.0    
    ##  [91] tibble_3.1.4           stringi_1.7.4          highr_0.9             
    ##  [94] RSpectra_0.16-0        forcats_0.5.1          lattice_0.20-44       
    ##  [97] Matrix_1.3-4           vctrs_0.3.8            pillar_1.6.2          
    ## [100] lifecycle_1.0.0        BiocManager_1.30.16    data.table_1.14.0     
    ## [103] bitops_1.0-7           R6_2.5.1               latticeExtra_0.6-29   
    ## [106] C50_0.1.5              rio_0.5.27             codetools_0.2-18      
    ## [109] MASS_7.3-54            gtools_3.9.2           assertthat_0.2.1      
    ## [112] openssl_1.4.5          withr_2.4.2            GenomeInfoDbData_1.2.6
    ## [115] hms_1.1.0              grid_4.1.1             rpart_4.1-15          
    ## [118] tidyr_1.1.3            rmarkdown_2.10         inum_1.0-4            
    ## [121] carData_3.0-4          partykit_1.2-15        WGCNA_1.70-3          
    ## [124] base64enc_0.1-3
