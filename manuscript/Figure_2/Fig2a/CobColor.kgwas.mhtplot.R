######################################################
###############chromosome lengths##################
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
#############draw mahtplot#############
######################################################
   gwas.data <- read.table("kmer.gwas.CC.p0.0001_perfect_align_to_B73_v4.txt",sep="\t",header=TRUE)
   p <- gwas.data$p
   mht.data.pro <- gwas.data[,c(2,3,4)]
   dimnames(mht.data.pro)[[2]] <- c("chr","pos","p")
# modify data structure.
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
   pngfile.names <- paste("Fig2a.CC_kmer_GWAS_p0.0001_B73_v4.png",sep="")
   title.name <- "Cob Color kmer GWAS"
   label.name<-c(paste("chr.",1:10,sep=""))
   color.array<-rep(c("red","blue"),5)
  
#open a png file for plotting.
   png(pngfile.names,res=200,height=1100,width=1500)
   y <- -log10(mht.data[,3])
   x <- mht.data[,2]
   plot(x,y,type="p",cex=0,xlab="Maize Chromosomes",ylab="",xlim=c(0,pos[10]),
        ylim=c(4,(max(y)+1)),xaxs="i",yaxs="i",xaxt="n",cex.lab=1.6, cex.axis=1.6)       
#use a 'for' cycle to draw points with differnt colors.
   for(k in 1:10){
       x.zc <- mht.data[mht.data[,1]==k,2]
       y.zc <- -log10(mht.data[mht.data[,1]==k,3])
       points(x.zc,y.zc,col=color.array[k],pch=19,cex=0.5)
   }      
#modify and embellish the mhtplot.
   axis(1,at=pos.chr,labels=paste(1:10,sep=""),cex.lab=1.6, cex.axis=1.6)
   dev.off()

