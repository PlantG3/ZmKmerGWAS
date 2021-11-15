library("readr")
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
outdir <- args[2]

### load KOC
d <- read_delim(file = input, delim = "\t")

### load normalization factor:
normfactors <- read.delim("3o-m282.k31.norm.factors.txt", comment.char = "#")

### remove low depth lines
low.depth.lines <- normfactors[normfactors$Remove == "yes", 1]
d <- d[, !colnames(d) %in% low.depth.lines]

normfactors <- normfactors[normfactors$Remove == "no", ]
nf <- normfactors$Normfactor
names(nf) <- normfactors$Geno
nfv <- as.numeric(nf[colnames(d)[-1]]) # normalization factor values

### correlation between RC and normalization count values
allcor <- apply(d[, -1], 1, cor, y = nfv)

### remove low variable kmers
d <- d[allcor <= 0.8, ]

### kmer normalization:
d2 <- round(t(t(d[, -1]) / nfv), 0) ### normalization
d2 <- cbind(d[, 1], d2) ### add Kmer

### output
out <- gsub(".*\\/", "", input)
out <- gsub("counts", "normKOC", out)
out <- paste0(outdir, "/", out)
write.table(d2, out, row.names=F, quote=F, sep="\t")
cat(input, "has been processed.\n")

