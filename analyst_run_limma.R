#Navin Ramanan - project 3
#Role - data analyst
#4/11/2022


#BiocManager::install("edgeR")

libs <- c("tidyverse", "ggVennDiagram", "BiocManager",
          "DESeq2", "edgeR", "limma")
# if you don't have a package installed, use BiocManager::install() or 
# install.packages(), as previously discussed.
#install.packages("statmod")
library(statmod)
library(limma)
library(ggplot2)

for (package in libs) {
  suppressPackageStartupMessages(require(package, 
                                         quietly = T, 
                                         character.only = T))
  require(package, character.only = T)
}

# sample info dataframe with array_id and chemical columns
read_samples_file <- function(sample_csv) {
  samples <- read.csv(sample_csv,as.is=TRUE)
  return(samples)
}

read_rma_matrix <- function(rma_matrix) {
  # the full RMA normalized matrix for all experiments
  rma <- read.table(rma_matrix,
                    sep='\t',
                    as.is=TRUE,
                    header=TRUE,
                    row.names=1,
  )
  return(rma)
}
  

# sample info dataframe with array_id and chemical columns
generate_limma_results <- function(samples, rma, sample_subset, vehicle) {

  sample2 <- filter(samples, samples$chemical %in% sample_subset)
  sample1 <- filter(sample2, sample2$vehicle == vehicle[1])
  # subset the full expression matrix to just those in sample1 comparison
  rma.subset <- rma[paste0('X',sample1$array_id)]
  # construct a design matrix modeling treatment  vs control for use by limma
  design <- model.matrix(~factor(sample1$chemical, levels=sample_subset))
  colnames(design) <- sample_subset 
  colnames(design)[1]<- 'Intercept'
  v = voom(rma.subset, design)
  fit = lmFit(v, method='robust', maxit=500)
  fit <- eBayes(fit)
  t <- topTable(fit, coef=2, n=nrow(rma.subset))
  # write out the results to file
  limma_trmt_file <- paste0('limma_',colnames(design)[2], '.csv')
  write.csv(t,limma_trmt_file)
  
  return(limma_trmt_file)
}


#### load and filter ####
read_limma_results <- function(filename) {
  limma_results <- (read.csv(filename, header = TRUE, sep = ','))
  limma_padj_results <- limma_results[order(limma_results$adj.P.Val, decreasing = TRUE), ]
  limma_padj_file <- paste0('padj_', filename)
  write.csv(limma_padj_results, limma_padj_file)
  message(paste0("Limma DE analysis result file sorted by adjusted p-value: ", 
                  limma_padj_file))
  
  return(limma_padj_results)
}


FC_hist_plot <- function(limma_result, title) {
  #histogram
  log2fc_plot <- ggplot(limma_result, aes(logFC))+
                 geom_histogram(color="black", fill="lightblue", bins=100) +
                 theme_bw() +
                 xlab("logFC") + 
                 ylab("count") +
                 labs(title = title)

  return(log2fc_plot)
}


DE_scat_plot <- function(limma_result, title) {
  #Scatter Plot
  scat_plot <- ggplot(limma_result, aes(x=logFC, y=P.Value)) + 
      geom_point(size=1.0) + 
      theme_bw() +
      labs(title = title) +
      xlab("LogFC") + ylab("P Value") 
  return(scat_plot)
  
}

#====================Concordance=======================
# convert refseq to gene
refseq_to_gene <- function(x, refseq_mapping){
  x <- na.omit(x)
  geneL <- c()
  for(refseq in x){
    gene <- refseq_mapping %>%
      dplyr::filter(REFSEQ == refseq) %>% # uses refseq ID
      dplyr::select(SYMBOL)
    gene <- unique(gene[,1])
    if(length(gene)==1){
      geneL <- c(gene, geneL)
    } else{
      gene <- NA
      geneL <- c(gene, geneL)
    }
  }
  geneL <- na.omit(geneL)
  geneL <- unique(geneL)
  return(geneL)
}

#Find the genes that have the same DE directionality.
agreement <- function(mic_data, rna_data){
  inter <- intersect(mic_data$REFSEQ, rna_data$X)
  agreed <- c()
  for(element in inter){
    sL <- sign(mic_data$logFC[which(mic_data$REFSEQ == element)])
    sL <- mean(sL)
    sD <- sign(rna_data$log2FoldChange[which(rna_data$X == element)])
    sD <- mean(sD)
    if(sL == sD){
      agreed <- c(element, agreed)
    }
  }
  return(agreed)
}

#Map probeid to refseq
probeid_mapping <- function(limma_data, affy_map) {

  names(limma_data)[1] <- "PROBEID"
  limma <- merge(affy_map, limma_data, by = "PROBEID", all.y = TRUE)
  limma <- subset (limma, select = -c(SYMBOL, GENENAME))
  limma <- limma[!duplicated(limma[c('PROBEID')]), ]
  limma <- limma[!is.na(limma$REFSEQ),]
  return(limma)
}

# Caiculate concordance score
calculate_concordance_score <- function(n0, n1, n2, N){
  nx <- ((n0*N)-(n1*n2))/(n0+N-n1-n2)
  x <-  (2 * nx) / (n1 + n2)
  y <-  (2 * n0) / (n1 + n2)
  return(y)
}  

#Plot Concordance vs DEs
concordance_vs_de_plot <- function(plot_data, x_label) {
  
  trmt_names <- c("3ME", "CLO", "CHL")
  de_plot <-ggplot(plot_data, aes(x=x_values, y=y_values, color=trmt_names)) +
            geom_point() +
            #add linear trend line
            geom_smooth(method=lm, se=FALSE, col='black', linetype="dashed", size=0.5) +
            theme_bw() + 
            scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
            xlab(x_label) + ylab('Concordance Percantage') +
            ggtitle(paste0("Concordance vs number of ", x_label, " for 3ME, CLO and CHL"))
  
  return(de_plot)
}


#===========================================================
FC_vs_Concordance_Plot <- function(data, m0, m1, m2, M){
  print("FC_ vs_Concordance_Plot")
  return(NULL)
  f.x <- function(n0,n1=m1,n2=m2,N=M) {
    x <- (n0*N-n1*n2)/(n0+N-n1-n2)
    if(x <= n0) {
      x
    } else {
      NA
    }
  }
  
  do.plot <- function(n1) {
    n0 <- seq(1,n1)
    N <- 20000
    plot(n0,
         vapply(n0,function(x) f.x(x,n1,100,N),0),
         type='l', lwd=2,
         xlab='Observed overlap', ylab='True intersection',
         col='red',
         theme_bw(),
         xlim=c(0,n1),
         ylim=c(0,n2)
    )
    lines(n0,
          vapply(n0,function(x) f.x(x,n1,1000,N),0),
          type='l', lwd=2,
          col='blue'
    )
    lines(n0,
          vapply(n0,function(x) f.x(x,n1,5000,N),0),
          type='l', lwd=2,
          col='green'
    )
    lines(n0,
          vapply(n0,function(x) f.x(x,n1,10000,N),0),
          type='l', lwd=2,
          col='grey'
    )
    lines(n0,n0,col="grey",lty=2)
    lines(n0,rep(0,length(n0)),col="grey",lty=2)
    legend(1,n1,legend=c(
      "Identity/Zero line",
      paste0("n1=",n1,", n2=100"),
      paste0("n1=",n1,", n2=1000"),
      paste0("n1=",n1,", n2=5000"),
      paste0("n1=",n1,", n2=10000")
    ),
    col=c("grey","red","blue","green","grey"),
    lty=c(2,1,1,1,1)
    )
  }
  print("n-0")
  print(x)
  fig3 <- do.plot(100)
  return(fig3)
}

##################################################################
