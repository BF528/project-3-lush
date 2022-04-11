library(reshape2)
library(ggplot2)
library(DESeq2)

#function to make count matrix from 9 featureCounts output files
DESeqDataSetFromFeatureCounts <- function (sampleTable, directory = ".", design, ignoreRank = FALSE, ...) 
{
  if (missing(design)) 
    stop("design is missing")
  l <- lapply(as.character(sampleTable[, 2]), function(fn) read.table(file.path(directory, fn), skip=2))
  if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1)))) 
    stop("Gene IDs (first column) differ between files.")
  tbl <- sapply(l, function(a) a$V7)
  colnames(tbl) <- sampleTable[, 1]
  rownames(tbl) <- l[[1]]$V1
  rownames(sampleTable) <- sampleTable[, 1]
  dds <- DESeqDataSetFromMatrix(countData = tbl, colData = sampleTable[, 
                                                                       -(1:2), drop = FALSE], design = design, ignoreRank, ...)
  return(dds)
}

#                                           
file_path <- '~/downloads/counts/'

#read in a manually assembled experimental design table 
sampleTable <- read.csv(paste0(file_path,'sample_table.txt'))

#create and extract count matrix
compiled <- DESeqDataSetFromFeatureCounts(sampleTable,directory = file_path,design = ~ condition)

counts <- compiled@assays@data$counts

#plot counts distribution across samples
boxplot <- ggplot(data = melt(as.data.frame(log10(counts)),id.vars=NULL),aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  scale_x_discrete(name = "Samples") +
  scale_y_continuous(name = "Counts") +
  labs(title="Distribution of Counts Across Samples") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#extract control sample counts from full count matrix
controls <- c('SRR1178000','SRR1178005','SRR1178007','SRR1177973','SRR1178016','SRR1178019')
control_table <- read.csv(paste0(file_path,'control_counts.csv'))
genes <- control_table[,1]

cords <- c()
for (i in controls) {
  cords <- c(cords,which(i==colnames(control_table)))
}

controls_filtered <- control_table[,cords]
rownames(controls_filtered) <- genes

#order rows in control and treated sample count matrices identically so they can be merged
#only merge if row names control==treated
ordered_controls <- controls_filtered[ order(row.names(controls_filtered)), ]
ordered_counts <- counts[ order(row.names(counts)), ]
if (setequal(rownames(ordered_counts),rownames(ordered_controls))==TRUE) {
  combined <- cbind(ordered_counts,ordered_controls)
  
}

#filter out genes with <10 counts across all samples
combined <- combined[rowSums(combined[])>10,]

conditions <- read.csv(paste0(file_path,'toxgroup_1_rna_info.csv'))

#get sample names for each treatment group and their associated controls
coldata <- conditions[,c(1,2,4)]

for_corn <- coldata[c(which(coldata$mode_of_action=='Control')),]

car <- coldata[c(which(coldata$mode_of_action=='CAR/PXR')),]
car <- car$Run

cyto <- coldata[c(which(coldata$mode_of_action=='Cytotoxic')),]
cyto <- cyto$Run

corn_controls <- for_corn[c(which(for_corn$vehicle=='CORN_OIL_100_%')),]
corn_controls <- corn_controls$Run
  
cmc <- coldata[c(which(coldata$vehicle=='CMC_.5_%')),]
cmc <- cmc$Run

#merge lists of treatment group sample names and their associated controls
ahr_de <- cmc
cyto_de <- c(cyto,corn_controls)
car_de <- c(car,corn_controls)

#subset treated sample and associated control counts from full count matrix
ahr_mat <- combined[,c(ahr_de)]
cyto_mat <- combined[,c(cyto_de)]
car_mat <- combined[,c(car_de)]

#prepare coldata variable for DESeq2
cords <- c()
for (i in ahr_de) {
  cords <- c(cords,which(i==coldata$Run))
}
ahr_coldata <- coldata[c(cords),]

cords <- c()
for (i in cyto_de) {
  cords <- c(cords,which(i==coldata$Run))
}
cyto_coldata <- coldata[c(cords),]

cords <- c()
for (i in car_de) {
  cords <- c(cords,which(i==coldata$Run))
}
car_coldata <- coldata[c(cords),]

#run DESeq2 for each treatment group
ahr_dds <- DESeqDataSetFromMatrix(countData = ahr_mat,
                              colData = ahr_coldata,
                              design = ~ mode_of_action)

cyto_dds <- DESeqDataSetFromMatrix(countData = cyto_mat,
                              colData = cyto_coldata,
                              design = ~ mode_of_action)

car_dds <- DESeqDataSetFromMatrix(countData = car_mat,
                              colData = car_coldata,
                              design = ~ mode_of_action)

#get differential expression results sorted by adjusted p-value and top 10 DE genes for treatment group
ahr_dds <- DESeq(ahr_dds)
ahr_res <- results(ahr_dds)
ahr_res <- ahr_res[ order(ahr_res$padj), ]
top_ahr <- rownames(ahr_res)[1:10]
write.csv(ahr_res,paste0(file_path,'ahr.csv'))

cyto_dds <- DESeq(cyto_dds)
cyto_res <- results(cyto_dds)
cyto_res <- cyto_res[ order(cyto_res$padj), ]
top_cyto <- rownames(cyto_res)[1:10]
write.csv(cyto_res,paste0(file_path,'cyto.csv'))

car_dds <- DESeq(car_dds)
car_res <- results(car_dds)
car_res <- car_res[ order(car_res$padj), ]
top_car <- rownames(car_res)[1:10]
write.csv(car_res,paste0(file_path,'car.csv'))

modes <- rep(c('AhR,3-METHYLCHOLANTHRENE','Cytotoxic','CAR/PXR'),each=10)
top_genes <- c(top_ahr,top_cyto,top_car)
top_tab <- cbind(top_genes,modes)
colnames(top_tab) <- c('Gene','mode_of_action')
write.csv(top_tab,paste0(file_path,'top_genes.csv'))

#get number of significantly differentially expressed genes per treatment group
low_padj <- c(length(which(ahr_res$padj<0.05)),
              length(which(cyto_res$padj<0.05)),
              length(which(car_res$padj<0.05)))

padj_tab <- cbind(unique(modes),low_padj)
colnames(padj_tab) <- c('mode_of_action','padj<0.05')
write.csv(padj_tab,paste0(file_path,'padjs.csv'))

#plot histograms of fold change per treatment group
par(mfrow=c(2,2))
AhR_log2FC <- ahr_res$log2FoldChange
hist_a <- hist(AhR_log2FC)
Cytotoxic_log2FC <- cyto_res$log2FoldChange
hist_cy <- hist(Cytotoxic_log2FC)
CARPXR_log2FC <- car_res$log2FoldChange
hist_ca <- hist(CARPXR_log2FC)

#plot fold change versus nominal p-value per treatment group
par(mfrow=c(2,2))
AhR_pvalue <- ahr_res$pvalue
plot(AhR_log2FC,AhR_pvalue, pch = 19, col = "blue")
Cytotoxic_pvalue <- cyto_res$pvalue
plot(Cytotoxic_log2FC,Cytotoxic_pvalue, pch = 19, col = "blue")
CARPXR_pvalue <- car_res$pvalue
plot(CARPXR_log2FC,CARPXR_pvalue, pch = 19, col = "blue")



