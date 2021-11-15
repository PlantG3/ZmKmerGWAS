x <- read.delim("DGAT.KOC.NAM.mean.txt")
a <- 0
b <- 0
pdf("Fig4e.oil_kmers_on_B73.pdf",width = 6.5,height = 5.5)
par(mar = c(5,5,2,5))
plot(a,b,xlim=c(1,length(x[,1])),ylim=c(min(x[,2:3]),max(x[,2:3])),xlab = "NAM parental lines",ylab = "Mean of KOCs",cex=0, cex.lab=2, cex.axis=2)
points(c(1:length(x[,1])),x$chr6_peak,col="black",pch=19,cex=1)
points(c(1:length(x[,1])),x$X469ins,col="red",pch=24,cex=1.5,lwd=2)
legend("topleft",
       legend=c("chr6 oil k-mers","F469ins k-mers"),ncol=1,
       lty=c(0,0),lwd=1,pch=c(19,24),col=c("black","red"),cex=1.5,x.intersp = 0.1)
dev.off()
