######################################################
###############The variables will be used##################
######################################################
chr.length <- c(301354135,237068873,232140174,241473504,217872852,
                169174353,176764762,175793759,156750706,150189435)
#use a 'for' cycle to produce several useful vectors.
for(i.chr in 1:10){
    if(i.chr==1){
        pos.elment1 <- sum(chr.length[1:i.chr])
        pos.elment2 <- pos.elment1/2
        pos <- pos.elment1
        pos.chr <- pos.elment2
    }else{
        pos.elment2 <- pos.elment1+chr.length[i.chr]/2
        pos.elment1 <- sum(chr.length[1:i.chr])
        pos <- c(pos,pos.elment1)
        pos.chr <- c(pos.chr,pos.elment2)}
}

######################################################
#############Using 'for' cycles to draw mahtplot#############
######################################################
   gwas.data <- read.table("Oil_NAM_RIL_QTL_results.txt",sep="\t",header=TRUE)
   #install.packages("randomcoloR")
   library(randomcoloR)
   f <- as.character(unique(gwas.data$RIL))
   n <- length(f)
   palette<-distinctColorPalette(n)
   cc <- cbind(f,palette)
   colnames(cc) <- c("RIL","color")
   gwas.data.color <- merge(gwas.data,cc,by="RIL")
   write.table(gwas.data.color,file="Oil_NAM_RIL_QTL_results_add_color.txt",sep="\t",quote=FALSE,row.names = FALSE)
   gwas.data <- gwas.data.color[order(gwas.data.color$chr,gwas.data.color$pos_v4),]
   p <- gwas.data$lod
   mht.data.pro <- gwas.data[,c(3,4,7)]
   dimnames(mht.data.pro)[[2]] <- c("chr","pos","p")
#use a second 'for' cycle to modify data structure.
   for(j in 2:10){
       if(j==2){
           chr.data <- mht.data.pro[mht.data.pro[,1]==j,]
           col2.data <- chr.data[,2]+pos[j-1]
           chr.data[,2] <-  col2.data
           mht.data <- rbind(mht.data.pro[mht.data.pro[,1]==1,],chr.data)
       }else{
           chr.data <- mht.data.pro[mht.data.pro[,1]==j,]
           col2.data <- chr.data[,2]+pos[j-1]
           chr.data[,2] <-  col2.data
           mht.data <- rbind(mht.data,chr.data)
            }
   }
#set some variables for mhtplot.
   color.array <- as.character(gwas.data$color)
  
#open a png file for plotting.
   pdf("Fig4c.oil_NAM_QTL_peaks.pdf",height = 5,width = 8)
   par(mar=c(5.1,5,4.1,1.2))
   y <- mht.data[,3]
   x <- mht.data[,2]
   plot(x,y,type="p",cex=0,xlab="Chromosome",ylab="LOD",xlim=c(0,pos[10]),
        ylim=c(0,(max(y)+5)),xaxs="i",yaxs="i",xaxt="n",cex.lab=2,cex.axis=2)
   
#use a 'for' cycle to draw points with differnt colors.
   for(k in 1:length(x)){
       points(x[k],y[k],col=color.array[k],pch=24,cex=1.5,lwd=2)
   }      
#modify and embellish the mhtplot.
   axis(1,at=pos.chr,labels=c(1:10),cex.lab=2,cex.axis=2)
   for (j in 1:9)
   {
     lines(x=c(sum(chr.length[1:j]),sum(chr.length[1:j])),y=c(0,max(y)+5),type="l",col="gray",lwd=2)
   }
   dev.off()
   pdf("figure_legend_new.pdf",height = 5,width = 8)
   par(mar=c(5.1,4.3,4.1,1.9))
   y <- mht.data[,3]
   x <- mht.data[,2]
   plot(x,y,type="p",cex=0,xlab="Chromosome",ylab="LOD",xlim=c(0,pos[10]),
        ylim=c(0,(max(y)+1)),xaxs="i",yaxs="i",xaxt="n",cex.lab=1.5,cex.axis=1.5)
   legend("topright",
          legend=f,ncol=3,
          lty=c(rep.int(0,23)),lwd=1,pch=c(rep.int(24,23)),col=color.array,cex=1,x.intersp = 0.1)
   dev.off()
