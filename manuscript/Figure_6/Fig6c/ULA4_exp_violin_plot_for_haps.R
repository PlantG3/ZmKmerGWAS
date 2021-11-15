library(ggplot2)
library(patchwork)
x <- read.delim("ULA4_NAM_exp_norm.txt")
head(x)
pdf("Fig5c.ULA4_exp_NAM_4_tissues_voilin_plot_for_haps.pdf",height = 8,width = 8)

### Ear
i <- 4
hap1 <- x[which(x$haplotype=="hap1"),i]
hap2 <- x[which(x$haplotype=="hap2"),i]
hap3 <- x[which(x$haplotype=="hap3"),i]
t.test(hap1,hap2)
t.test(hap1,hap3)
t.test(hap2,hap3)
data <- data.frame(name=c(rep("Hap1",length(hap1)),rep("Hap2",length(hap2)),rep("Hap3",length(hap3))),value=c(hap1,hap2,hap3))
p <- ggplot(data, aes(x=name, y=value,fill=name)) +
  geom_violin(trim=FALSE)+
  labs(title=colnames(x)[i],x="Haplotypes", y = "Gene Expression")+
  geom_boxplot(width=0.1,fill=c("white"))+
  theme_classic()
p1 <- p+scale_fill_manual(values=c("red", "blue","green"))

### Root
i <- 5
hap1 <- x[which(x$haplotype=="hap1"),i]
hap2 <- x[which(x$haplotype=="hap2"),i]
hap3 <- x[which(x$haplotype=="hap3"),i]
t.test(hap1,hap2)
t.test(hap1,hap3)
t.test(hap2,hap3)
data <- data.frame(name=c(rep("Hap1",length(hap1)),rep("Hap2",length(hap2)),rep("Hap3",length(hap3))),value=c(hap1,hap2,hap3))
p <- ggplot(data, aes(x=name, y=value,fill=name)) +
  geom_violin(trim=FALSE)+
  labs(title=colnames(x)[i],x="Haplotypes", y = "Gene Expression")+
  geom_boxplot(width=0.1,fill=c("white"))+
  theme_classic()
p2 <- p+scale_fill_manual(values=c("red", "blue","green"))

### Shoot
i <- 6
hap1 <- x[which(x$haplotype=="hap1"),i]
hap2 <- x[which(x$haplotype=="hap2"),i]
hap3 <- x[which(x$haplotype=="hap3"),i]
t.test(hap1,hap2)
t.test(hap1,hap3)
t.test(hap2,hap3)
data <- data.frame(name=c(rep("Hap1",length(hap1)),rep("Hap2",length(hap2)),rep("Hap3",length(hap3))),value=c(hap1,hap2,hap3))
p <- ggplot(data, aes(x=name, y=value,fill=name)) +
  geom_violin(trim=FALSE)+
  labs(title=colnames(x)[i],x="Haplotypes", y = "Gene Expression")+
  geom_boxplot(width=0.1,fill=c("white"))+
  theme_classic()
p3 <- p+scale_fill_manual(values=c("red", "blue","green"))

### Tassel
i <- 7
hap1 <- x[which(x$haplotype=="hap1"),i]
hap2 <- x[which(x$haplotype=="hap2"),i]
hap3 <- x[which(x$haplotype=="hap3"),i]
t.test(hap1,hap2)
t.test(hap1,hap3)
t.test(hap2,hap3)
data <- data.frame(name=c(rep("Hap1",length(hap1)),rep("Hap2",length(hap2)),rep("Hap3",length(hap3))),value=c(hap1,hap2,hap3))
p <- ggplot(data, aes(x=name, y=value,fill=name)) +
  geom_violin(trim=FALSE)+
  labs(title=colnames(x)[i],x="Haplotypes", y = "Gene Expression")+
  geom_boxplot(width=0.1,fill=c("white"))+
  theme_classic()
p4 <- p+scale_fill_manual(values=c("red", "blue","green"))

p1+p2+p3+p4

dev.off()