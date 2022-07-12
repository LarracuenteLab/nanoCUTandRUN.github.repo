library(ggplot2)
library(tidyverse)
library(dplyr)

pdf(paste0("Figure3.pdf"),width=12, height=10, paper="a4r")

par(mfrow=c(2,2))

for(repeatID in c("HETA","TAHRE", "TART-B1", "TART-A") ){
  
#for hiphop
protein="hiphop"

control <- read.table(paste("WT_1_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
colnames(control) <- c("position", "control_RPM")

treatment_rep1 <- read.table(paste(protein,"_1_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
colnames(treatment_rep1) <- c("position", "treatment_RPM1")

treatment_rep2 <- read.table(paste(protein,"_2_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
colnames(treatment_rep2) <- c("position", "treatment_RPM2")

HipHop_Rep1=(treatment_rep1[,2]+1)/(control[,2]+1)
HipHop_Rep2=(treatment_rep2[,2]+1)/(control[,2]+1)
hiphop=cbind(control[,1], HipHop_Rep1, HipHop_Rep2)

Mean=(hiphop[,2]+hiphop[,3])/2
hiphop=cbind(hiphop, Mean)
colnames(hiphop)=c("position", "Rep1","Rep2", "hiphop Mean")



#for hoap
protein="hoap"

control <- read.table(paste(protein,"_WG_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
colnames(control) <- c("position", "control_RPM")

treatment_rep1 <- read.table(paste(protein,"_1_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
colnames(treatment_rep1) <- c("position", "treatment_RPM1")

treatment_rep2 <- read.table(paste(protein,"_2_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
colnames(treatment_rep2) <- c("position", "treatment_RPM2")


Hoap_Rep1=(treatment_rep1[,2]+1)/(control[,2]+1)
Hoap_Rep2=(treatment_rep2[,2]+1)/(control[,2]+1)
Hoap=cbind(control[,1], Hoap_Rep1, Hoap_Rep2)

Mean=(Hoap[,2]+Hoap[,3])/2
Hoap=cbind(Hoap, Mean)
colnames(Hoap)=c("position", "Rep1","Rep2", "Hoap Mean")

plot(hiphop[,1],hiphop[,4], type="l", col="darkgreen", main=repeatID, xlab="Position (bp)", ylab="Target/Control")
lines(Hoap[,1],Hoap[,4], type="l", col="darkblue")
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", c("hiphop", "HOAP"), lty=1, lwd=3, col=c("darkgreen", "darkblue"), bty="n", ncol=2, cex=1.5)

dev.off()



pdf(paste0("FigureS6.pdf"),width=12, height=10, paper="a4r")
par(mfrow=c(2,2))
for(repeatID in c("Heta-1_D","Heta-2","Heta-3", "Heta-5") ){
  

  #for hiphop
  protein="hiphop"
  
  control <- read.table(paste("WT_1_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
  colnames(control) <- c("position", "control_RPM")
  
  treatment_rep1 <- read.table(paste(protein,"_1_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
  colnames(treatment_rep1) <- c("position", "treatment_RPM1")
  
  treatment_rep2 <- read.table(paste(protein,"_2_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
  colnames(treatment_rep2) <- c("position", "treatment_RPM2")
  
  
  HipHop_Rep1=(treatment_rep1[,2]+1)/(control[,2]+1)
  HipHop_Rep2=(treatment_rep2[,2]+1)/(control[,2]+1)
  hiphop=cbind(control[,1], HipHop_Rep1, HipHop_Rep2)
  
  Mean=(hiphop[,2]+hiphop[,3])/2
  hiphop=cbind(hiphop, Mean)
  colnames(hiphop)=c("position", "Rep1","Rep2", "hiphop Mean")
  
  
  
  #for hoap
  protein="hoap"
  
  control <- read.table(paste(protein,"_WG_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
  colnames(control) <- c("position", "control_RPM")
  
  treatment_rep1 <- read.table(paste(protein,"_1_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
  colnames(treatment_rep1) <- c("position", "treatment_RPM1")
  
  treatment_rep2 <- read.table(paste(protein,"_2_",repeatID,"_dimer_vs_monomer.distrib", sep=""), header=TRUE)
  colnames(treatment_rep2) <- c("position", "treatment_RPM2")
  
  
  Hoap_Rep1=(treatment_rep1[,2]+1)/(control[,2]+1)
  Hoap_Rep2=(treatment_rep2[,2]+1)/(control[,2]+1)
  Hoap=cbind(control[,1], Hoap_Rep1, Hoap_Rep2)
  
  Mean=(Hoap[,2]+Hoap[,3])/2
  Hoap=cbind(Hoap, Mean)
  colnames(Hoap)=c("position", "Rep1","Rep2", "Hoap Mean")
  
  plot(hiphop[,1],hiphop[,4], type="l", col="darkgreen", main=repeatID, xlab="Position (bp)", ylab="Target/Control")
  lines(Hoap[,1],Hoap[,4], type="l", col="darkblue")
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", c("hiphop", "HOAP"), lty=1, lwd=3, col=c("darkgreen", "darkblue"), bty="n", ncol=2, cex=1.5)

dev.off()

