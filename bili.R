setwd('D:/R/home/rate histogram')
library("stringr")
library("do")
temp = list.files(pattern="*.txt")
for (names in temp)
{
name = strsplit(names,'.',1)[[1]][1]
print(name)
tables = read.table(str_c("D:/R/home/rate histogram/",name,".txt"),sep="\t",header=T)
tables$Status = Replace(tables$Status,"1", "Dead")
tables$Status = Replace(tables$Status,"0", "Alive")
library(tidyverse)
tables <- count(tables,class,Status)

pdf(str_c("D:/R/home/rate histogram/",name,".pdf"),width=length(unique(tables$class))*2,height = 8)

library("ggplot2")
row_num <- length(row(tables))/length(tables)
temp_pct <-matrix(tables$n, nrow = row_num)


for(i in 1:(row_num/2))
{
  temp_sum = temp_pct[i*2-1,1] + temp_pct[i*2,1]
  temp_pct[i*2-1,1] <- temp_pct[i*2-1,1]/temp_sum
  temp_pct[i*2,1] <- temp_pct[i*2,1]/temp_sum
}
print(tables)
p <- ggplot(tables,aes(fill=Status, y= temp_pct, x = class,fill=variable)) + 
  geom_col()+
  geom_text(aes(label = scales::percent(temp_pct)))+
  scale_y_continuous(labels = scales::percent) + 
  labs(y = "Rate",x=name) +
  theme_gray()
print(p)
dev.off()
}