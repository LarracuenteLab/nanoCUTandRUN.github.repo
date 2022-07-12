pdf("FigureS6.pdf", paper="a4", width = 12, height = 12)

par(mfrow=c(2,1))
data=read.csv("hoap.CenIsland_MergeHETA.common.peak", sep="\t")

row.names(data)=data[,1]

par(mar=c(6, 5, 4, 1))
plot(data[1:20,2],col="red", xaxt="n", ylab= "Normalized RPM count (IP/Control)", xlab="", main="HOAP", ylim=c(0,120), pch=16)
points(data[1:20,3],col="red", pch=16)
axis(1,1:20,labels=FALSE)

text(1:20,rep(-10,20),labels=rownames(data)[1:20],srt = 45,xpd=NA,adj=c(1,1), cex=0.6)



data=read.csv("hiphop.CenIsland_MergeHETA.common.peak", sep="\t")

row.names(data)=data[,1]

par(mar=c(6, 5, 4, 1))

plot(data[1:20,3], xaxt="n", ylab= "Normalized RPM count (IP/Control)",col="red", xlab="", main="HipHop", ylim=c(0,130), pch=16)
points(data[1:20,2],col="red", pch=16)
axis(1,1:20,labels=FALSE)

text(1:20,rep(-10,20),labels=rownames(data)[1:20],srt = 45,xpd=NA,adj=c(1,1), cex=0.6)

dev.off()


##### Merging table
data_HipHop=read.csv("hiphop.common.peak", sep="\t")
colnames(data_HipHop)=c("Repeat", "HipHop Rep2 RPM count (IP/Control)", "HipHop Rep1 RPM count (IP/Control)")

data_Hoap=read.csv("hoap.common.peak", sep="\t")
colnames(data_Hoap)=c("Repeat", "HOAP Rep2 RPM count (IP/Control)", "HOAP Rep1 RPM count (IP/Control)")

data=merge(data_HipHop, data_Hoap, by="Repeat", all=T)
data=data[order(data$`HipHop Rep2 RPM count (IP/Control)` , decreasing = T),]
write.table(data, "Merge.table.csv", quote=F, row.names = F, col.names = T, sep=",")
