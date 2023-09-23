
library(dplyr)
library(corrplot)

setwd("~/cor_tri")

my_data <- read.csv('~/machine_learning/data/tcga/survival_out.csv',header=T,row.names=1,check.names=F)


res <- cor(my_data)
#method = "square"
pdf("1.pdf",width=length(colnames(my_data))/3,height=length(colnames(my_data))/3)
corrplot(res, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 90,col=myColor)
dev.off()