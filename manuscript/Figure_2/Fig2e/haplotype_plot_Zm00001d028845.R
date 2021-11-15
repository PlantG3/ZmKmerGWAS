x <- read.delim("Zm00001d028845_SNPs.txt")
y <- x[,-1:-5]
c <- c()
for (i in 1:length(y[,1]))
{
  maj <- which(y[i,]==as.character(x$major[i]))
  min <- which(y[i,]==as.character(x$minor[i]))
  het <- which(y[i,]=="H")
  mis <- which(y[i,]=="N")
  a <- c()
  a[maj] <- "red"
  a[min] <- "blue"
  a[het] <- "green"
  a[mis] <- "gray"
  c <- rbind(c,a)
}
a <- 0
b <- 0
i <- 1
j <- 1
pdf("Fig2d.haplotype_plot_Zm00001d028845.pdf")
plot(a,b,xlim=c(1,length(y[,1])),ylim = c(0,length(y[1,])),xlab="Marker",ylab="Line",cex=0,main="Zm00001d028845")
for (i in 1:length(y[,1]))
{
  for (j in 1:length(y[1,]))
  {
    points(i,j,col=as.character(c[i,j]),pch=19,cex=0.6)
  }
}
dev.off()