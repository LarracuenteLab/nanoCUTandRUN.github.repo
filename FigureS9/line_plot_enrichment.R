library(ggplot2)
library(tidyverse)
library(dplyr)

#repeatID=c("HETA", "Heta-1_D", "Heta-2","Heta-3", "Heta-5", "HETRP", "TAHRE", "TART-A", "TART-B1", "TART-C_NTPR")
#for(i in c(1:10)){

#c("HETA","Heta-1_D","Heta-2","Heta-3","Heta-5","TAHRE", "TART-B1", "TART-A", "TART-C_NTPR")

pdf(paste0("FigureS7.pdf"),width=15, height=15, paper="a4")

par(mfrow=c(4,2))

for(repeatID in c("HETA","Heta-1_D","Heta-2","Heta-3","Heta-5","TAHRE", "TART-B1", "TART-A") ){
  
#repeatID="TART-C_NTPR"
 
hiphop <- read.table(paste("hiphop_WG_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
colnames(hiphop) <- c("position", "RPM")

hoap <- read.table(paste("hoap_WG_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
colnames(hoap) <- c("position", "RPM")

WT <- read.table(paste("WT_WG_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
colnames(WT) <- c("position", "RPM")

aa=max(hiphop[,2], hoap[,2], WT[,2])
plot(hiphop[,1],hiphop[,2], type="l", col="darkgreen", main=repeatID, xlab="Position (bp)", ylab="Whole genome coverage (RPM)", ylim=c(0,aa))
lines(hoap[,1],hoap[,2], type="l", col="darkblue")
#lines(WT[,1],WT[,2], type="l", col="darkred")

}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("hiphop", "HOAP"), lty=1,lwd=3, col=c("darkgreen", "darkblue"), bty="n", ncol=3, cex=1.2)

dev.off()

