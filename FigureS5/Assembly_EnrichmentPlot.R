library(karyoploteR)
library(regioneR)
library(GenomicRanges) 
library(rtracklayer) 
library(IRanges) 
library(devtools)
library(stringr)

# Building a Genome
mygenome= read.table("dmel_scaffold2_plus0310.chrom.sizes",comment.char = "")
mygenome2 = data.frame(mygenome[,1],rep(1,dim(mygenome)[1]),mygenome[,2])
custom.genome <- toGRanges(mygenome2)


# Building a Cytoband
cytobands = read.csv("dmel_scaffold2_plus0310_2.Jockey3.cytoband.mod", sep="\t", header=T)
custom.cytobands <- toGRanges(cytobands)

# Peak
peaks_Hiphop2<-read.table("hiphop-2min_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks_Hiphop2  <- toGRanges(peaks_Hiphop2)

peaks_Hiphop15<-read.table("hiphop-15min_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks_Hiphop15  <- toGRanges(peaks_Hiphop15)

peaks_Hoap2<-read.table("hoap-2min_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks_Hoap2  <- toGRanges(peaks_Hoap2)

peaks_Hoap15<-read.table("hoap-15min_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks_Hoap15  <- toGRanges(peaks_Hoap15)


# Colors Cytoband

color_legend=read.csv("name-color-telomere.csv", sep=";", header=F)
vv=as.character(color_legend$V2)
names(vv)=color_legend$V1

chr = "4_2"
detail.region <- toGRanges(data.frame(chr, 1000000,1317735))

pdf("plot/4_2.pdf")

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

#AA="visible.region"
AA=500
BB=500

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp, zoom=detail.region)
#kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)

kp<- kpPlotBigWig(kp, data="hiphop-15min.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop-15min.mapped.sorted.q30.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)

#BB <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0 , ymax=BB,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop-15min", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop15, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hiphop-2min.mapped.sorted.RPM.bw", ymax=AA ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop-2min.mapped.sorted.q30.RPM.bw", ymax=AA ,r0=0.66, r1=0.79, col="black",border=NA)
#BB <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop-2min", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop2, data.panel = 1, r0=0.64,r1=0.65,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap-15min.mapped.sorted.RPM.bw", ymax=AA, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap-15min.mapped.sorted.q30.RPM.bw", ymax=AA, r0=0.49, r1=0.62, col="black",border=NA)
#BB <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0 , ymax=BB, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP-15min", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap15, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap-2min.mapped.sorted.RPM.bw", ymax=AA, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap-2min.mapped.sorted.q30.RPM.bw", ymax=AA, r0=0.32, r1=0.45, col="black",border=NA)
#BB <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0 , ymax=BB, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP-2min", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap2, data.panel = 1, r0=0.30,r1=0.31,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="WT-15min.mapped.sorted.RPM.bw", ymax=BB, r0=0.17, r1=0.32, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT-15min.mapped.sorted.q30.RPM.bw", ymax=BB, r0=0.17, r1=0.32, col="black",border=NA)
#BB <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0 , ymax=BB, r0=0.15, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT-15min",  r0=0.15, r1=0.28, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT-2min.mapped.sorted.RPM.bw", ymax=BB, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT-2min.mapped.sorted.q30.RPM.bw", ymax=BB, r0=0, r1=0.13, col="black",border=NA)
#BB <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0 , ymax=BB, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT-2min",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()



