##################################################################################################
# input
##################################################################################################
datapath <- "."
genome <- "A188"
trait <- "4f3o_ULA"
chr.size <- read.delim("A188Ref1.length", header=F, stringsAsFactors=F)
colnames(chr.size) <- c("Chr", "Size")
row.names(chr.size) <- chr.size$Chr
chromosomes <- 1:10
chr.size <- chr.size[chr.size$Chr %in% chromosomes, ]
chr.size <- chr.size[as.character(chromosomes), ]
xlab.text <- "chromosomes"
ylab.text <- substitute(paste("-log10(", italic(p), ")"))
background_cols <- c("grey90", "grey80")
cols <- c("olivedrab4", "lightsalmon3")

##################################################################################################
# data
##################################################################################################
samfiles <- dir(path=datapath, pattern=paste0("mapping.", genome, ".sam"), full.names=T)
neg_samfiles <- grep("negative", samfiles, value=T)
pos_samfiles <- grep("positive", samfiles, value=T)

possam <- read.delim(pos_samfiles, header=F, stringsAsFactors=F)
possam <- possam[possam[,3] %in% chromosomes, ]
possam <- possam[, c(1, 3, 4)]
colnames(possam) <- c("Seq", "Chr", "Pos")
possam$Pval <- gsub(".*_", "", possam[,1])
possam$neglogP <- -log10(as.numeric(possam$Pval))

negsam <- read.delim(neg_samfiles, header=F, stringsAsFactors=F)
negsam <- negsam[negsam[,3] %in% chromosomes, ]
negsam <- negsam[, c(1, 3, 4)]
colnames(negsam) <- c("Seq", "Chr", "Pos")
negsam$Pval <- gsub(".*_", "", negsam[,1])
negsam$neglogP <- -log10(as.numeric(negsam$Pval))

##################################################################################################
# plot
##################################################################################################
pdffile <- paste0(trait, ".kmerGWAS.", genome, ".sigplot.pdf")
pdf(pdffile, width=5, height=3.5)

accum <- 0
all.chr <- NULL
centers <- NULL
gap <- sum(chr.size$Size)/80
ymax <- max(c(possam$neglogP, negsam$neglogP))
ymin <- min(c(possam$neglogP, negsam$neglogP))
yrange <- c(ymin, ymax)
yextend <- (ymax - ymin) / 50

### initiate plotting
par(mar=c(3.5, 3.5, 3, 0.5), mgp=c(2, 0.5, 0))
plot(NULL, NULL, ylim = yrange,
     xlim = c(0, gap * (nrow(chr.size)-1) + sum(chr.size$Size)),
     xaxt = "n", yaxt = "n",
     xlab = xlab.text, ylab = ylab.text,
     main = paste(trait, genome, sep="-"), cex.lab = 1,
     bty = "n")

axis(side = 2, cex.axis = 1, lwd = 1)
box(lwd = 1)

all.accum <- NULL
if (length(cols) < nrow(chr.size)) {
  cols <- rep(cols, ceiling(nrow(chr.size) / length(cols)))
}

accum <- 0
all.col <- NULL
all.chr <- NULL
centers <- NULL
for (i in 1:(nrow(chr.size))) {
  all.accum <- c(all.accum, accum)
  pre.accum <- accum
  chr <- chr.size[i, "Chr"]
  len <- chr.size[i, "Size"]
  
  # plot shades
  polygon(x=c(accum, accum, accum+len, accum+len), border=NA,
          y=c(yrange[1]-yextend, yrange[2]+yextend,
              yrange[2]+yextend, yrange[1]-yextend),
          col=background_cols[i %% 2 + 1])
  
  # positive points
  if (sum(possam$Chr %in% chr) > 0) {
    pos <- as.numeric(possam[possam$Chr %in% chr, "Pos"])
    prob <- possam[possam$Chr %in% chr, "neglogP"]
    points(accum+pos, prob, lwd=0.6, pch=21, cex=0.6, col="black", bg=cols[1])
  }
  
  # negative points
  if (sum(negsam$Chr %in% chr) > 0) {
    pos <- as.numeric(negsam[negsam$Chr %in% chr, "Pos"])
    prob <- negsam[negsam$Chr %in% chr, "neglogP"]
    points(accum+pos, prob, lwd=0.6, pch=21, cex=0.6, col="black", bg=cols[2])
  }
  
  accum <- accum + len + gap
  center.point <- (pre.accum + accum - gap)/2
  #all.col <- c(all.col, plot.col)
  all.chr <- c(all.chr, chr)
  centers <- c(centers, center.point)
}

# axis
axis(side=1, at=centers, labels=all.chr, tick=F,
     cex.axis=1, lwd=1)

dev.off()
