####loading data####
geno <- read.table("1o_ZeaGBSv2.7_B73xB97.Z001E.NAM.RILs.geno.txt",header = T,sep = "\t",stringsAsFactors = F)
map <- read.table("1o_ZeaGBSv2.7_B73xB97.Z001E.genetic.map",header = T)
phe <- read.csv("282_NAM_RIL_NIROil_mean.csv",header = T)

####convert data format####
geno.map <- cbind(geno[,1],map[,4],geno[,2:length(geno[1,])])
colnames(geno.map)[1:2] <- c("rs","gmap")
chr <- strsplit(as.character(geno.map$rs),"_")
chr <- unlist(lapply(chr, function(x) x[1]))
chr <- gsub("S","",chr)
geno.map.chr <- cbind(geno.map[,1],chr,geno.map[,2:length(geno.map[1,])])
geno.phe.name <- cbind(colnames(geno)[-1],c(1:length(colnames(geno)[-1])))
colnames(geno.phe.name) <- c("Line","flag")
geno.phe <- merge(phe,geno.phe.name,by="Line")
flag <- as.numeric(as.vector(geno.phe$flag))
geno.map.chr.phe <- geno.map.chr[,c(1:3,flag+3)]
geno.map.chr.phe.t <- t(geno.map.chr.phe)
row.names(geno.map.chr.phe.t)[1:3] <- c("Genotype","","")
write.table(geno.map.chr.phe.t,file="ZeaGBSv2.7_B73xB97.Z001E.format.txt",sep="\t",quote=FALSE,col.names = FALSE)
geno.phe <- geno.phe[,2:1]
colnames(geno.phe) <- c("UpperLeafAngle","Genotype")
write.csv(geno.phe,file="ZeaGBSv2.7_B73xB97.Z001E.phe.csv",row.names = FALSE)
csv <- read.table("ZeaGBSv2.7_B73xB97.Z001E.format.txt",sep="\t",header = TRUE)
write.csv(csv,file="ZeaGBSv2.7_B73xB97.Z001E.format.csv",row.names = FALSE)

############################QTL start############################
#install.packages("qtl",dependencies = TRUE)
############################QTL start############################
library(qtl)
## read the table ##
phe <- read.csv("ZeaGBSv2.7_B73xB97.Z001E.phe.csv",header = T,sep = ",")
geno <- read.csv("ZeaGBSv2.7_B73xB97.Z001E.format.csv",header = T,sep = ",")
sug <- read.cross("csvs", dir = ".", genfile = "ZeaGBSv2.7_B73xB97.Z001E.format.csv", phefile = "ZeaGBSv2.7_B73xB97.Z001E.phe.csv",genotypes = c("A","H","B"))
sug <- jittermap(sug)
## calculate prob ##
prosug <- calc.genoprob(sug)
out.em <- scanone(prosug,model = "normal",method = "em")
out.em_df <- as.data.frame(out.em)
out.em_v2 <- cbind(row.names(out.em_df),out.em_df)
colnames(out.em_v2) <- c("marker","chr","pos_g","lod")
out.em_v3 <- merge(out.em_v2,map,by="marker")
out.em_v3 <- out.em_v3[,c(1,2,6,3,4)]
colnames(out.em_v3)[2] <- "chr"
write.table(out.em_v3,file=paste0("ZeaGBSv2.7_B73xB97.Z001E.lod.txt"),sep="\t",quote=FALSE,row.names = FALSE)
## permutation test to ensure the LOD cutoff, n=1000##
operm <- scanone(prosug, method="em", n.perm=1000)
summary(operm, alpha=c(0.05,0.1,0.2))
qtl <- summary(out.em, perms=operm, alpha=0.05, pvalues=TRUE)
qtl <- cbind(row.names(qtl),qtl)
colnames(qtl)[1] <- "rs"
write.table(qtl,file="ZeaGBSv2.7_B73xB97.Z001E.QTL.txt",sep="\t",quote=FALSE,row.names = FALSE)
pdf(file="ZeaGBSv2.7_B73xB97.Z001E.QTL.pdf",height=5,width = 8)
plot(out.em)
dev.off()
## fitqtl to get the pve
qc <- as.numeric(qtl$chr)
qp <- round(as.numeric(qtl$pos))
qtl <- makeqtl(prosug, chr=qc, pos=qp, what="prob")
out.fqi <- fitqtl(prosug, qtl=qtl, method="hk",formula = y~Q1+Q2)## if we have two QTL
sum.qtl <- out.fqi$result.drop
sum.qtl <- cbind(row.names(sum.qtl),sum.qtl)
colnames(sum.qtl)[1] <- "rs"
write.table(sum.qtl,file="ZeaGBSv2.7_B73xB97.Z001E.QTL.PVE.txt",sep="\t",quote=FALSE,row.names = FALSE)