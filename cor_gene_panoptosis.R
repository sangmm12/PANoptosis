setwd("D:/R/home/cor")
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplotify)


tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

out_put_list <- c()
for (tumor_id in tumor_list)
{
print(tumor_id)
dat_input <- read.csv("D:/R/home/cor/every/panoptosis.csv",check.names = F,header=T)

dat_input <- dat_input[grep(tumor_id,dat_input$CODE),]
dat_input <- dat_input[grep("Tumor",dat_input$Group),]
dat_input <- dat_input[,-1:-2]

dat_input <- dat_input[!duplicated(dat_input$SampleName),]

rownames(dat_input) <- make.unique(sub('\n','',dat_input$SampleName))

dat_input <- dat_input[,-(length(colnames(dat_input)))]

colscluster = length(colnames(dat_input))/3


data_all_tumor <- read.csv(paste("D:/R/home/cor/every/CIBER.csv",sep=''),check.names = F)

dat_one_tumor <- data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]

if (length(colnames(dat_one_tumor))==0)
{
    next
}

rownames(dat_one_tumor) <- make.unique(sub('\n','',dat_one_tumor$SampleName))

dat_one_tumor <- dat_one_tumor[,!duplicated(colnames(dat_one_tumor))]

dat_one_tumor <- dat_one_tumor[,3:length(colnames(dat_one_tumor))]



#CIBER
if (tumor_id == "CHOL"||tumor_id=="PCPG")
{
  dat_one_tumor <- select(dat_one_tumor,-c('T_cells_gamma_delta'))
}
if (tumor_id == "LUAD"||tumor_id=="PCPG")
{
  dat_one_tumor <- select(dat_one_tumor,-c('T_cells_CD4_naive'))
}
if (tumor_id=="TGCT")
{
  dat_one_tumor <- select(dat_one_tumor,-c("Eosinophils"))
}


one_tumor_sample <- unlist(rownames(dat_one_tumor))

all_name <- names(which(table(c(rownames(dat_input),one_tumor_sample))==2))

dat_gene <- dat_one_tumor[match(all_name,rownames(dat_one_tumor)),]

dat_im <- dat_input[match(all_name,rownames(dat_input)),]
print(dat_im)
library(psych)
data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")

data.r <- data.corr$r
data.p <- data.corr$p

write.csv(data.p,file=paste("out_p.csv",sep=''),quote=F)
library(pheatmap)


getSig <- function(dc) {
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}

sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)


paletteLength <- 1000
myColor <- colorRampPalette(c("green", "white", "red"))(paletteLength)



test <- data.r
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))



pdf(paste("pancancer/panoptosis/",tumor_id,".pdf",sep=''),width = length(colnames(dat_input))/1.5,height = length(colnames(data_all_tumor))/3)

pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,cluster_cols=F,cellwidth = 20,cellheight = 20, display_numbers=sig.mat,fontsize=length(colnames(dat_input)))

dev.off()
}