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
peaks_Hiphop1<-read.table("hiphop_1.unique_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks_Hiphop1<-toGRanges(peaks_Hiphop1)

peaks_Hoap1<-read.table("hoap_1.unique_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks_Hoap1<- toGRanges(peaks_Hoap1)


# Colors Cytoband

color_legend=read.csv("name-color-telomere.csv", sep=";", header=F)
vv=as.character(color_legend$V2)
names(vv)=color_legend$V1

################# ################# ################# 
################# Tel X ################# 
################# ################# ################# 

pdf("plot_1Rep/telX.pdf")

chr = "X_1"
detail.region <- toGRanges(data.frame(chr, 1,123800))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

AA=800
BB=40
CC=800
DD=40
EE=800
FF=40

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)

kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)


kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)


dev.off()


################# ################# ################# 
################# Tel 2L ################# 
################# ################# ################# 

pdf("plot_1Rep/tel2L.pdf")

chr = "2L_1"
detail.region <- toGRanges(data.frame(chr, 1,35000))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

AA=800
BB=40
CC=800
DD=40
EE=800
FF=40

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=5000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)


kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)


kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Tel 2R ################# 
################# ################# ################# 

pdf("plot_1Rep/tel2R.pdf")

chr = "2R_21"
detail.region <- toGRanges(data.frame(chr, 21579400,21691270))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10


AA=800
BB=40
CC=800
DD=40
EE=800
FF=40

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)


kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)


kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Tel 3L_1 ################# 
################# ################# ################# 

pdf("plot_1Rep/tel3L.pdf")

chr = "3L_1"
detail.region <- toGRanges(data.frame(chr, 1,350000))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

AA=50
BB=20
CC=350
DD=20
EE=200
FF=20

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)


kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)


kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Tel 3R ################# 
################# ################# ################# 

pdf("plot_1Rep/tel3R.pdf")

chr = "3R_28"
detail.region <- toGRanges(data.frame(chr, 27800000,27994903))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

AA=100
BB=20
CC=300
DD=20
EE=100
FF=20

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)

kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)


kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Tel 4 ################# 
################# ################# ################# 

pdf("plot_1Rep/4_2.pdf")

chr = "4_2"
detail.region <- toGRanges(data.frame(chr, 1000000,1317735))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

AA=700
BB=40
CC=700
DD=40
EE=700
FF=20

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)


kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)


kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)

dev.off()

################# ################# ################# 
################# Legend ################# 
################# ################# ################# 

pdf("plot_1Rep/legend_telomere.pdf")
plot.new()
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("center",legend=color_legend[,1], pch=15,col=vv,cex=0.8, ncol=3)

dev.off()

################# ################# ################# 
################# Cen 4 ################# 
################# ################# ################# 
# Colors Cytoband

color_legend=read.csv("name-color-centromere.csv", sep=";", header=F)
vv=as.character(color_legend$V2)
names(vv)=color_legend$V1


pdf("plot_1Rep/Cen4.pdf")

chr = "Contig119"
#detail.region <- toGRanges(data.frame(chr, 1000000,1317735))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

AA=20
BB=20
CC=20
DD=20
EE=20
FF=20

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)


kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)


kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Cen X ################# 
################# ################# ################# 

pdf("plot_1Rep/CenX.pdf")

chr = "Contig79"
#detail.region <- toGRanges(data.frame(chr, 1000000,1317735))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

AA=80
BB=100
CC=80
DD=100
EE=80
FF=100

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)


kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)


kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Cen 3 ################# 
################# ################# ################# 

pdf("plot_1Rep/Cen3.pdf")

chr = "3R_5"
detail.region <- toGRanges(data.frame(chr, 17300,75000))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

AA=10
BB=20
CC=10
DD=20
EE=10
FF=20

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)


kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)


kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Cen 2 ################# 
################# ################# ################# 

pdf("plot_1Rep/Cen2.pdf")

chr = "tig00057289"
#detail.region <- toGRanges(data.frame(chr, 17300,75000))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

AA=10
BB=10
CC=10
DD=10
EE=10
FF=10

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)


kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)


kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)

dev.off()



################# ################# ################# 
################# Cen Y ################# 
################# ################# ################# 

pdf("plot_1Rep/CenY.pdf")

chr = "Y_Contig26"
#detail.region <- toGRanges(data.frame(chr, 17300,75000))

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

AA=10
BB=10
CC=10
DD=10
EE=10
FF=10

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)

kp<- kpPlotBigWig(kp, data="hiphop_1.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_1.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.64, r1=0.89, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.64, r1=0.89, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop", r0=0.64, r1=0.89, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="hoap_1.mapped.sorted.RPM.bw", ymax=CC, r0=0.32, r1=0.57, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_1.unique.mapped.sorted.RPM.bw", ymax=CC,r0=0.32, r1=0.57, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.32, r1=0.57, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP", r0=0.32, r1=0.57, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="orange",border=NA)

kp<- kpPlotBigWig(kp, data="WT_1.mapped.sorted.RPM.bw", ymax=EE, r0=0, r1=0.25, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_1.unique.mapped.sorted.RPM.bw", ymax=EE,  r0=0, r1=0.25, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0, r1=0.25, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT",  r0=0, r1=0.25, cex=0.6, label.margin = 0.035)

dev.off()

