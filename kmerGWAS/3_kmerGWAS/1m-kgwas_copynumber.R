# input arguments:
args <- commandArgs(trailingOnly=T)
kset <- as.numeric(args[1])
meanDepth <- as.numeric(args[2])
indir <- args[3]
trait <- args[4]
pc <- args[5]
outdir <- args[6]
max.pval <- as.numeric(as.character(args[7]))

library("readr")

####### convert KOCs to copy number #############
count2cn <- function(count, cmean) {
	cn <- count
	nonzero <- cn[count>0]
	# for nonzero counts, set to >=1
	# squeeze copy number due to count variation
	nonzero_cn <- pmax(1, round(nonzero / cmean / 2))
	cn[count>0] <- nonzero_cn
	cn
}
#################################################


#################################################
lmgwas <- function(geno, mean_count, pheno, pc1, pc2, pc3) {
# function for k-mer GWAS
# number to copy number transformation
	pheno <- as.numeric(pheno)
	geno <- as.numeric(geno)
	geno <- count2cn(geno, mean_count)
	lmfull <- lm(pheno ~ geno + pc1 + pc2 + pc3)
	lmred <- lm(pheno ~ pc1 + pc2 + pc3)
	pval <- anova(lmred, lmfull)[2, 6]
	coeff <- lmfull$coefficients[2]
	c(coeff, pval)
}
#################################################

# data
phe <- read.csv(trait)
phe.lines <- phe[,1]
km.all <- system(paste("ls -1", indir), intern=T)
km <- km.all[kset]
cat(km, "\n")
kmer <- read_delim(paste0(indir, "/", km), delim="\t",
	col_types=cols(.default="d", Kmer=col_character()))
	kmerall <- kmer
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
out <- apply(kmer[,-1], 1, lmgwas, mean_count=meanDepth, pheno=phe[,2], pc1=pc1, pc2=pc2, pc3=pc3)

selected <- (!is.na(out[2, ]) & out[2, ]<=max.pval)  # only selected low-pvalue tests for output

coeff <- out[1, ]
pvals <- out[2, ]
coeff[!is.na(coeff)] <- format(coeff[!is.na(coeff)], digits=3, scientific=T)
pvals[!is.na(pvals)] <- format(pvals[!is.na(pvals)], digits=3, scientific=T)

# output
out <- cbind(kmer[,1], coeff, pvals)
trait <- colnames(phe)[2]
colnames(out) <- c("Kmer", paste(trait, c("coeff", "pvals"), sep="_") )
write.table(out,file=paste0(outdir, "/", km, ".", trait, ".kgwas.txt"), sep="\t", quote=F,row.names=F)

# KOC output
kocout <- cbind(out[selected, ], kmerall[selected, -1])
write.table(kocout,file=paste0(outdir, "/", km, ".", trait, ".kgwas.KOC.txt"), sep="\t", quote=F,row.names=F)

