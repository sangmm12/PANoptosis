rm(list=ls())
file_dir = "D:/R/home/classify_pancancer"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
need_list = c('Lasso')
for (file_name in need_list)
{
tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

final_tumor_list = c()

my_data = read.csv('yzx.csv',header=T,check.names=F)

other_data = read.csv(paste(file_name,".csv",sep=''),header=T,check.names=F)

p_value_csv <- data.frame()

for(name in tumor_list)
{

dat <- data.frame(check.names = F)

other_file <- other_data[grep(name,other_data$CODE),]

if (length(rownames(other_file))==0)
{
    next
}

other_file <- other_file[!duplicated(other_file$SampleName),]

rownames(other_file) = gsub('\n','',other_file[,length(other_file)])

other_file <- other_file[,2:3]

other_file <- other_file[colnames(other_file)[1]]

exp_file <- my_data[grep(name,my_data$CODE),]

exp_file <- exp_file[grep("Tumor",exp_file$Group),]

exp_file <- exp_file[,-1:-2]

exp_file <- exp_file[!duplicated(exp_file$SampleName),]

rownames(exp_file) = gsub('\n','',exp_file[,length(exp_file)])

exp_file <- exp_file[,1:length(exp_file)-1]

gene_list <- colnames(exp_file)

temp_length <- length(gene_list)

all_name <- names(which(table(c(rownames(other_file),rownames(exp_file)))==2))

if (length(all_name)==0)
{
    next
}

for(gene_name in gene_list)
{
    for(i in all_name)
    {
        dat <- rbind(dat,c(gene_name,other_file[match(i,rownames(other_file)),],exp_file[c(gene_name)][match(i,rownames(exp_file)),]))
    }

}

colnames(dat) <- c("Gene","Group","value")

dat[,3] = as.numeric(dat[,3])


xx <-compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")

p_value <- as.matrix(xx$p)

p_value[is.na(p_value)] <- 1



final_tumor_list <- append(final_tumor_list,name)
print(name)



pdf(paste(file_name,"/",name,".pdf",sep=''),width=temp_length,height = 8)

p <- ggboxplot(dat, x = "Gene", y = "value",
               palette = "simpsons", 
               color = "Group",x.text.angle=60,x.text.size=20)

p <- p + xlab("")+ylab("Expression Value")
p <- p + theme(axis.text = element_text(size = 15),axis.title=element_text(size=30))
print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))

dev.off()
#rm(list=c('dat','other_file','exp_file','all_name','xx','p_value','p','temp_length'))
}

}
