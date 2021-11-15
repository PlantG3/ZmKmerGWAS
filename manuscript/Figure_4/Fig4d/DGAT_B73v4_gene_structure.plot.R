# setting
geno <- "B73v4"
setwd(".")
colors <- c("olivedrab4", "lightsalmon3")
par(mgp=c(2, 0.6, 0), mar=c(3.5, 3.5, 3, 0.5))
pdfout <- "Fig4d.DGAT_B73v4_gene_strucutre.pdf"
# data
dataPath <- "B73v4_DGAT"
faslen <- read.delim(paste0(dataPath, "/B73v4_DGAT.2.Zm00001d036982.length"), header=F)
maxlen <- faslen[1,2]
beds <- dir(path=paste0(dataPath, "/B73v4_DGAT.4.pos.adjusted.gtf.bed"), pattern="bed$", full.names=T)
nbeds <- length(beds)

negative_map_sam <- "B73v4_DGAT_k.negative.mapping.sam"
positive_map_sam <- "B73v4_DGAT_k.positive.mapping.sam"

# plot
pdf(pdfout, width=6, height=5)
main_info <- paste0("oil k-mers on the DGAT1-2 (", geno, ")")
plot(NULL, NULL, xlim=c(1, maxlen), ylim=c(0, nbeds+1),
     xlab="position (bp)", ylab="", yaxt="n", cex.lab=1.2, cex.axis=1.2,
     main=main_info)
base <- 0.4
exonheight <- 0.4
gstart <- NULL
gend <- NULL

# gene structure
for (bed in beds) {
  lines(x=c(1, maxlen), y=c(base, base), col="gray80", lwd=2)
  beddata <- read.delim(bed, header=F)
  if (is.null(gstart)) {
    gstart <- min(beddata[,2])
  } else if (min(beddata[,2]) < gstart) {
    gstart <- min(beddata[,2])
  }
  
  if (is.null(gend)) {
    gend <- max(beddata[,3])
  } else if (max(beddata[,2]) > gend) {
    gend <- max(beddata[,3])
  }
  
  for (i in 1:nrow(beddata)) {
    color <- beddata[i, 7]
    rect(xleft=beddata[i,2], xright=beddata[i,3],
         ybottom=base-exonheight/2, ytop=base+exonheight/2,
         col=as.character(color), border=as.character(color))
  }
  
  # transcript information
  transcript <- gsub(".adjusted.bed", "", gsub(".*\\/", "", bed))
  text(0, base+exonheight, labels=transcript, pos=4)
  
  base <- base + 0.8
}

# indicate gene location
rect(xleft=gstart, xright=gend,
     ybottom=base, ytop=base+1.6,
     col="gray90", border=NA)

# assocKmers
neg_ak <- read.delim(negative_map_sam, header=F, flush=T)
pos_ak <- read.delim(positive_map_sam, header=F, flush=T)
neg_lowest_pval <- min(as.numeric(gsub(".*_", "", neg_ak[, 1])))
neg_highest_pval <- max(as.numeric(gsub(".*_", "", neg_ak[, 1])))
pos_lowest_pval <- min(as.numeric(gsub(".*_", "", pos_ak[, 1])))
pos_highest_pval <- max(as.numeric(gsub(".*_", "", pos_ak[, 1])))
lowest_pval <- min(neg_lowest_pval, pos_lowest_pval)
highest_pval <- min(neg_highest_pval, pos_highest_pval)
neglog10pval_bottom <- floor(-log10(highest_pval))
neglog10pval_top <- ceiling(-log10(lowest_pval))

base_bottom <- base - 0.2
base_top <- base + 1.6

# module to convert p-value to plot y position
pval2pos <- function(p, pbottom, ptop, pos_bottom, pos_top) {
  pos_bottom + (p - pbottom) * (pos_top - pos_bottom) / (ptop - pbottom)
}

# lines and axis of pvals
ypos_lines <- pval2pos(neglog10pval_bottom:neglog10pval_top, neglog10pval_bottom, neglog10pval_top, base_bottom, base_top)
abline(h=ypos_lines, col="gray80", lwd=0.8)
axis(2, at=ypos_lines, labels = neglog10pval_bottom:neglog10pval_top, las=2, cex.axis=1.2)
mtext(text=substitute(paste("-log10(", italic(p), ")")), side=2, line = 2, at=median(ypos_lines), cex=1.2)

# kmer plotting
for (i in 1:nrow(neg_ak)) {
  pval <- as.numeric(gsub(".*_", "", neg_ak[i, 1]))
  neglog10_pval <- -log10(pval)
  plot_ypos <- pval2pos(neglog10_pval, neglog10pval_bottom, neglog10pval_top, base_bottom, base_top)
  nmismatch <- as.numeric(gsub(".*\\:", "", neg_ak[i, 13]))
  #lines(rep(neg_ak[i, 4], 2), c(base+0.1, base+0.3), col=colors[2])
  shape <- 15
  if (nmismatch==0) {
    shape <- 19
  }
  points(neg_ak[i, 4], plot_ypos, pch=shape, col=colors[2])
}


for (i in 1:nrow(pos_ak)) {

  pval <- as.numeric(gsub(".*_", "", pos_ak[i, 1]))
  neglog10_pval <- -log10(pval)
  plot_ypos <- pval2pos(neglog10_pval, neglog10pval_bottom, neglog10pval_top, base_bottom, base_top)
  nmismatch <- as.numeric(gsub(".*\\:", "", pos_ak[i, 13]))
  #lines(rep(neg_ak[i, 4], 2), c(base+0.1, base+0.3), col=colors[2])
  shape <- 15
  if (nmismatch==0) {
    shape <- 19
  }
  points(pos_ak[i, 4], plot_ypos, pch=shape, col=colors[1])
}

dev.off()
