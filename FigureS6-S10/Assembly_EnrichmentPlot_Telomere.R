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
peaks_Hiphop1<-read.table("hiphop_2.unique_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks_Hiphop1<-toGRanges(peaks_Hiphop1)

peaks_Hoap1<-read.table("hoap_2.unique_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks_Hoap1<- toGRanges(peaks_Hoap1)


# Colors Cytoband

color_legend=read.csv("name-color-telomere.csv", sep=";", header=F)
vv=as.character(color_legend$V2)
names(vv)=color_legend$V1

################# ################# ################# 
################# Tel X ################# 
################# ################# ################# 

pdf("plot_Rep1_WG/telX.pdf")

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

kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Tel 2L ################# 
################# ################# ################# 

pdf("plot_Rep1_WG/tel2L.pdf")

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

kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Tel 2R ################# 
################# ################# ################# 

pdf("plot_Rep1_WG/tel2R.pdf")

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
kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Tel 3L_1 ################# 
################# ################# ################# 

pdf("plot_Rep1_WG/tel3L.pdf")

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

kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()

################# ################# ################# 
################# Tel 4 ################# 
################# ################# ################# 

pdf("plot_Rep1_WG/4_2.pdf")

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

kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()



################# ################# ################# 
################# Tel 3R ################# 
################# ################# ################# 

pdf("plot_Rep1_WG/tel3R.pdf")

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

kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()



################# ################# ################# 
################# Cen 4 ################# 
################# ################# ################# 
# Colors Cytoband

color_legend=read.csv("name-color-centromere.csv", sep=";", header=F)
vv=as.character(color_legend$V2)
names(vv)=color_legend$V1


pdf("plot_Rep1_WG/Cen4.pdf")

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

kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Cen X ################# 
################# ################# ################# 

pdf("plot_Rep1_WG/CenX.pdf")

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

kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Cen 3 ################# 
################# ################# ################# 

pdf("plot_Rep1_WG/Cen3.pdf")

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
EE=20
FF=20

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=1, chromosomes=chr, plot.params = pp, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=2, cex=0.7, digits=3 )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=8)

kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()


################# ################# ################# 
################# Cen 2 ################# 
################# ################# ################# 

pdf("plot_Rep1_WG/Cen2.pdf")

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

kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()



################# ################# ################# 
################# Cen Y ################# 
################# ################# ################# 

pdf("plot_Rep1_WG/CenY.pdf")

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

kp<- kpPlotBigWig(kp, data="hiphop_2.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_2.unique.mapped.sorted.RPM.bw", ymax=AA, r0=0.83, r1=0.96, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=AA,r0=0.83, r1=0.96, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop Rep2", r0=0.83, r1=0.96, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hiphop1, data.panel = 1, r0=0.81,r1=0.82,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.66, r1=0.79, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hiphop_WG.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hiphop_WG.unique.mapped.sorted.RPM.bw", ymax=BB ,r0=0.66, r1=0.79, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=BB, r0=0.66, r1=0.79, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HipHop WG", r0=0.66, r1=0.79, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="hoap_2.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_2.unique.mapped.sorted.RPM.bw", ymax=CC, r0=0.49, r1=0.62, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=CC, r0=0.49, r1=0.62, numticks = 2, cex=0.6)
kpAddLabels(kp, labels = "HOAP Rep2", r0=0.49, r1=0.62, cex=0.6, label.margin = 0.035)
kpPlotRegions(kp, data=peaks_Hoap1, data.panel = 1, r0=0.47,r1=0.48,avoid.overlapping=FALSE,col="orange",border=NA)

kpDataBackground(kp, r0=0.32, r1=0.45, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="hoap_WG.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="hoap_WG.unique.mapped.sorted.RPM.bw", ymax=DD, r0=0.32, r1=0.45, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=DD, r0=0.32, r1=0.45, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "HOAP WG", r0=0.32, r1=0.45, cex=0.6, label.margin = 0.035)

kp<- kpPlotBigWig(kp, data="WT_2.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_2.unique.mapped.sorted.RPM.bw", ymax=EE, r0=0.17, r1=0.28, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=EE, r0=0.17, r1=0.28, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT Rep2",  r0=0.17, r1=0.28, cex=0.6, label.margin = 0.035)

kpDataBackground(kp, r0=0, r1=0.13, color = "#FFFFE0")
kp<- kpPlotBigWig(kp, data="WT_WG.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data="WT_WG.unique.mapped.sorted.RPM.bw", ymax=FF, r0=0, r1=0.13, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=FF, r0=0, r1=0.13, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "WT WG",  r0=0, r1=0.13, cex=0.6, label.margin = 0.035)

dev.off()

