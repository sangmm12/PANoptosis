setwd('~/gene_distribution')

library(RCircos)

data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, tracks.inside=10, tracks.outside=0 )

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size <- 0.7
RCircos.Reset.Plot.Parameters(rcircos.params)

out.file<-"RCircosDemoHumanGenome.pdf"
pdf(file=out.file,height=8,width=8,compress=TRUE)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
data(RCircos.Gene.Label.Data) 
name.col <- 4
data = read.table('ML.txt',sep='\t')
RCircos.Gene.Connector.Plot(data,
                            + 1, "in")

RCircos.Gene.Name.Plot(data,
                       + name.col,2, "in")
dev.off()