library(ggplot2)
x <- read.csv("Zm00001d028845_RC.csv")
head(x)
r <- x[,grep("R_",colnames(x))]
w <- x[,grep("W_",colnames(x))]
r_rc <- apply(r, 2, mean)
w_rc <- apply(w,2,mean)
data <- data.frame(name=c(rep("Red",length(r_rc)),rep("White",length(w_rc))),value=c(r_rc,w_rc))
pdf("Fig2d.violin_plot_KOC_in_Zm00001d028845.pdf")
p <- ggplot(data, aes(x=name, y=value,fill=name)) + 
  geom_violin(trim=FALSE)+
  labs(title="Zm00001d028845",x="Cob Color", y = "Mean of CoKs")+
  geom_boxplot(width=0.1,fill=c("white"))+
  theme_classic()
p+scale_fill_manual(values=c("red", "white"))
dev.off()