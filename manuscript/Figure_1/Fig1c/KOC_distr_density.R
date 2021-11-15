x <- read.delim("kmer_depth_sum_interval_100.txt")
head(x)
a <- x$kmer_sum_interval[-1]
b <- 10^(x$all_kmer_log10_freq_sum)[-1]
pdf("Fig1c.distribution_of_KOCs.pdf")
plot(a,b,type="l",lwd=2,col="black",ylim=c(0,2E8),xaxt="n",yaxt="n",xlab="Normalized k-mer depth",ylab="Number of k-mer")
#points(a,b,col="blue",pch=19,cex=0.5)
axis(1,at=c(seq(from=0,to=10000,by=1000)),labels = c(paste0(0:9),">10"),cex.axis=1.3)
axis(2,at=c(seq(from=0,to=2E8,by=5E7)),labels = c(0,50,100,150,200),cex.axis=1.3)
dev.off()