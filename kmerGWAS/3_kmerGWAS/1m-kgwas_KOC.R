# input arguments:
args <- commandArgs(trailingOnly=T)
kset <- as.numeric(args[1])
indir <- args[2]
trait <- args[3]
pc <- args[4]
outdir <- args[5]
max.pval <- as.numeric(as.character(args[6]))

#log=paste0(outdir, ".log")
#cat("Kmer GWAS", "\n", file=log)
#cat(paste("Trait input   :", trait), "\n", file=log, append=T)
#cat(paste("Kmer directory:", indir), "\n", file=log, append=T)
#cat(paste("PCA input     :", pc), "\n", file=log, append=T)
#cat(paste("Kmer set      :", kset), "\n", file=log, append=T)
#cat(paste("Maximum pvalue:", max.pval), "\n", file=log, append=T)
#cat(paste("outpur directory", outdir), "\n", file=log, append=T)

library("readr")

#################################################
lmgwas <- function(geno, pheno, pc1, pc2, pc3) {
# function for k-mer GWAS
	pheno <- as.numeric(pheno)
	geno <- as.numeric(geno)
	lmfull <- lm(pheno ~ geno + pc1 + pc2 + pc3)
	lmred <- lm(pheno ~ pc1 + pc2 + pc3)
	pval <- anova(lmred, lmfull)[2, 6]
	coeff <- lmfull$coefficients[2]
	c(coeff, pval)
}
#################################################

# data
#flag.txt  kfiles.txt  PC3.txt  trait.txt
phe <- read.delim(trait)
phe.lines <- phe[,1]
km.all <- system(paste("ls -1", indir), intern=T)
km <- km.all[kset]
cat(km, "\n")
kmer <- read_delim(paste0(indir, "/", km), delim="\t",
	col_types=cols(.default="d", Kmer=col_character()))
kmer.lines <- colnames(kmer)[-1]
popstr <- read.delim(pc)
popstr.lines <- popstr[, 1]

# common lines
common.lines <- intersect(intersect(phe.lines, kmer.lines), popstr.lines)
cat(length(common.lines), "common lines identified\n")

phe <- phe[phe[,1] %in% common.lines, ]
kmer <- kmer[, c(colnames(kmer)[1], common.lines)]
popstr <- popstr[popstr[,1] %in% common.lines, ]

pc1 <- popstr$PC1; pc2 <- popstr$PC2; pc3 <- popstr$PC3

# kmer GWAS
out <- apply(kmer[,-1], 1, lmgwas, pheno=phe[,2], pc1=pc1, pc2=pc2, pc3=pc3)

selected <- (out[2, ] <= max.pval)  # only selected low-pvalue tests for output

coeff <- format(out[1, ], digits=3, scientific=T)
pvals <- format(out[2, ], digits=3, scientific=T)

# output
out <- cbind(kmer[,1], coeff, pvals)
trait <- colnames(phe)[2]
colnames(out) <- c("Kmer", paste(trait, c("coeff", "pvals"), sep="_") )
write.table(out,file=paste0(outdir, "/", km, ".", trait, ".kgwas.txt"), sep="\t", quote=F,row.names=F)

# KOC output
kocout <- cbind(out[selected, ], kmer[selected, -1])
write.table(kocout,file=paste0(outdir, "/", km, ".", trait, ".kgwas.KOC.txt"), sep="\t", quote=F,row.names=F)

