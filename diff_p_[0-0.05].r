
library("pheatmap")
library("jsonlite")

setwd(dir = "D:/R/home/heatmap/p[0-0.05]")
temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)
myfiles = lapply(myfiles, na.omit)
file_nums = length(temp)
filename = sapply(strsplit(temp,"\\."),"[[",1)
for(i in filename)
{
df = read.csv(paste0(i,'.csv'),header=T,row.names=1)
df = replace(df,is.na(df),1)
df_temp = df

df = -log10(abs(df))
df[df_temp<0] = -df[df_temp<0]

pdf(paste0(i,'.pdf'),length(colnames(df))/2,length(rownames(df))/2)


paletteLength = 1000

getSig <- function(dc) {
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''}
  return(sc)
}

sig.mat <- matrix(sapply(as.matrix(df_temp), getSig), nrow=nrow(as.matrix(df_temp)))
str(sig.mat)

xx <- pheatmap(df,
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F,cluster_cols=F, cellwidth = 20,cellheight = 20,main="-log10(p)",display_numbers=sig.mat)
print(xx)
dev.off()
}
