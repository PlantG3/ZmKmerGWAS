library(ggplot2)
y <- read.delim("DGAT2_2.k31.m282.teo.landrance.normKOC.txt")
f <- as.character(unique(y$module))
i <- 1
x <- y[which(y$module==f[i]),c(-1:-7)]
head(x)
t <- x[,grep("teo",colnames(x))]
l <- x[,grep("lr",colnames(x))]
m <- x[,-c(1,2,grep("teo",colnames(x)),grep("lr",colnames(x)))]
t_rc <- apply(t,2, mean)
l_rc <- apply(l,2,mean)
m_rc <- apply(m,2,mean)
data <- data.frame(name=c(rep("Teosinte",length(t_rc)),rep("Inbred",length(m_rc))),value=c(t_rc,m_rc))
pdf(file="Fig4f.DGAT1_2_KOC_maize282_teo.pdf",height = 5,width = 6)
p <- ggplot(data, aes(x=name, y=value,fill=name)) + ylim(c(0,4)) +
  geom_violin(trim=FALSE)+
  labs(title=paste0(f[i]," module"),x="Group", y = "Mean of KOC")+
  geom_point(position = position_jitter(seed = 1, width = 0.2))+
  theme_classic()
p+scale_fill_manual(values=c("red", "blue"))
dev.off()
