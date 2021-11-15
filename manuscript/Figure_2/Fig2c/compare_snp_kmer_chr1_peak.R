s <- read.delim("snp.gwas.CC.p0.0001.txt",sep="\t",header=TRUE)
k <- read.delim("kmer.gwas.CC.p0.0001_perfect_align_to_B73_v4",sep="\t",header=TRUE)
sp <- s[which(s$Chr==1 & s$Pos_v4 >= 47691734 & s$Pos_v4 <= 48930655),]
head(sp)
kp <- k[which(k$chr==1 & k$pos >= 47691734 & k$pos <= 48930655),]
head(kp)
a <- 0
b <- 0
pdf("Fig2c.compare_snp_kmer_chr1_peak.pdf")
plot(a,b,xlim=c(47691734,48930655),ylim=c(min(-log10(c(sp$p,kp$p))),max(-log10(c(sp$p,kp$p)))),xlab="Position",ylab="p-value",main="Compare SNP and k-mer GWAS peak")
points(kp$pos,-log10(kp$p),col="blue",pch=19,cex=0.3)
points(sp$Pos_v4,-log10(sp$p),col="red",pch=19,cex=0.3)

legend("topleft", legend=c("SNP","k-mer"), bty = "n",
      pch=rep(19, 2), col=c("red","blue"),cex=1.5)
dev.off()