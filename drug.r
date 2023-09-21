rm(list=ls())
file_dir = "D:/R/home/classify_pancancer"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)

need_list = c('Lasso')
for (file_name in need_list)
{
tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

final_tumor_list = c()

other_data = read.csv(paste(file_name,".csv",sep=''),header=T,check.names=F)

p_value_csv <- matrix(nrow = 28,ncol = 8,  byrow = TRUE)
colnames(p_value_csv) = c('Bortezomib', 'Dactinomycin', 'Docetaxel', 'Daporinad', 'Sepantronium bromide', 'Vinblastine', 'Staurosporine','Dinaciclib')
rownames(p_value_csv) = tumor_list

for(name in tumor_list)
{

exp_file = read.csv(paste('drug/',name,'.csv',sep = ''),header=T,check.names=F,row.names = 1)

exp_file <- exp_file[,-length(exp_file)]

gene_list = colnames(exp_file)

temp_length <- length(gene_list)








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

for (p_name in colnames(exp_file))
{
    if (p_name%in%colnames(p_value_csv))
    {
        p_value_csv[which(tumor_list==name),c(p_name)] = p_value[which(colnames(exp_file)==p_name)]
    }
}

print(name)

pdf(paste(file_name,"/",name,".pdf",sep=''),width=temp_length,height = 8)
p <- ggboxplot(dat, x = "Gene", y = "value",
               palette = "ucscgb", 
               fill = "Group",x.text.angle=60,x.text.size=20)

p <- p + xlab("Drug")+ylab("sensitivity (IC50)")
p <- p + theme(axis.text = element_text(size = 15),axis.title=element_text(size=30))

write.table(compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova",label = "p.signif")$p.signif,file=paste(file_name,"/",name,".txt",sep=''),row.names=F,col.names=F,sep='\n',quote=F)
print(p + stat_compare_means(aes(group = Group),label = "p.signif",method = "anova"))


dev.off()

}

write.csv(p_value_csv,file=paste(file_name,"/p.csv",sep=''),quote=F)

}