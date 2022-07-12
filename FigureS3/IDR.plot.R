
Tel=c("4_2", "X_1", "3R_28", "3L_1", "2R_21", "2L_1")
Cen=c("Contig79", "Contig119", "3R_5", "tig00057289", "Y_Contig26")


data_HipHop = read.table("HipHop_TopCommonPeaks.bed")

Telomere = sum(summary(data_HipHop$V1)[Tel], na.rm=T)
Centromere = sum(summary(data_HipHop$V1)[Cen], na.rm=T)
Other=dim(data_HipHop)[1]-(Telomere+Centromere)


pdf("peak.localization.HipHop.pdf")
barplot(c(Telomere, Centromere, Other), names.arg=c("Telomere", "Centromere", "Other"), col=c("darkmagenta", "darkred", "bisque4"), main="HipHop", ylim=c(0,1800), ylab="peak number")
dev.off()


data_Hoap = read.table("Hoap_TopCommonPeaks.bed")

Telomere = sum(summary(data_Hoap$V1)[Tel], na.rm=T)
Centromere = sum(summary(data_Hoap$V1)[Cen], na.rm=T)
Other=dim(data_Hoap)[1]-(Telomere+Centromere)


pdf("peak.localization.HOAP.pdf")
barplot(c(Telomere, Centromere, Other), names.arg=c("Telomere", "Centromere", "Other"), col=c("darkmagenta", "darkred", "bisque4"), main="HOAP", ylim=c(0,250), ylab="peak number")
dev.off()


