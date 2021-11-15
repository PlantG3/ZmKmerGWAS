library(ggplot2)
x <- read.delim("ccd1_KOC.txt")
y <- x[,grep("Y_",colnames(x))]
w <- x[,grep("W_",colnames(x))]
a <- 0
b <- 0
ym <- apply(y,2,mean)
wm <- apply(w,2,mean)
data <- data.frame(name=c(rep("Yellow",length(ym)),rep("White",length(wm))),value=c(ym,wm))
pdf("Fig3g.ccd1_KOC_white_yellow.pdf")
p <- ggplot(data, aes(x=name, y=value,fill=name)) + 
  geom_violin(trim=FALSE)+
  labs(title="Zm00001d048373 (ccd1)",x="Kernel Color", y = "Mean depth of k-mers")+
  geom_boxplot(width=0.1,fill=c("white","yellow"))+
  theme_classic()
p+scale_fill_manual(values=c("white", "yellow"))
dev.off()