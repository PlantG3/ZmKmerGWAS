#######################################################################
# read KOCs of oil k-mers
koc <- read.delim("4h3o_oil.k31.m282.teo.landrance.normKOC.txt", check.names=F, stringsAsFactors=F)
head(koc)
lines <- colnames(koc)[-(1:3)]
landrace <- grep("lr_", lines, value=T)
teosinte <- grep("teo_", lines, value=T)
m282 <- lines[!grepl("lr_|teo_", lines) ]

#######################################################################
chisq <- function(x, maize_cols, teo_cols) {
# chisq test to test the null hypothesis
# no difference in k-mer presence in maize and teosinte populations
  maize <- as.numeric(as.character(x[maize_cols]))
  teo <- as.numeric(as.character(x[teo_cols]))
  nma <- sum(maize==0)
  nmp <- sum(maize>0)
  nta <- sum(teo==0)
  ntp <- sum(teo>0)
  #if ((nma+nmp)==0 | (nta+ntp)==0) {
  #  cat(x, "\n")
  #}
  
  chidata <- matrix(c(nma, nmp, nta, ntp), nrow=2, byrow=T)
  chipval <- chisq.test(chidata)[[3]]
  nmp_freq <- nmp / (nma + nmp)
  ntp_freq <- ntp / (nta + ntp)
  c(chipval, nmp_freq, ntp_freq)
}
#######################################################################
maize_cols <- which(colnames(koc) %in% m282)
teo_cols <- which(colnames(koc) %in% teosinte)
chiout <- apply(koc, 1, chisq, maize_cols=maize_cols, teo_cols=teo_cols)

kout <- koc[, 1:3]
kout$maize_teo_freq_Pval <- chiout[1, ]
kout$maize_freq <- chiout[2, ]
kout$teo_freq <- chiout[3, ]

#pval_cutoff <- 0.05 / nrow(kout) # Bonferroni control
padj <- p.adjust(kout$maize_teo_freq_Pval, method="BH")  ### FDR method
pval_cutoff <- max(kout$maize_teo_freq_Pval[which(padj==max(padj[padj<0.05]))])

nrow(kout)
sum(kout$maize_teo_freq_Pval < pval_cutoff, na.rm=T)

kout_sig <- kout[!is.na(kout$maize_teo_freq_Pval) & kout$maize_teo_freq_Pval < pval_cutoff, ]

# output selected kmers
write.table(kout_sig, "4h4o_historic.selection.kmers.txt", quote=F, row.names=F, sep="\t")

kout_sig <- read.delim("4h4o_historic.selection.kmers.txt")
###############################################################################################################
pdf("Fig4g.oil_selected_kmers.pdf", width=4.5, height=4.5)

par(mgp=c(2, 0.5, 0))
plot(kout_sig$maize_freq - kout_sig$teo_freq, kout_sig$NIROil_coeff,
     xlab="% maize - % teosinte (k-mer presence)", ylab="K-mer coeff to oil",
     cex=0.5, col="gray30", main="Evolutionarily selected oil k-mers")

dgat <- kout_sig[which(kout_sig$label_P39=="c6_112.06M"),]
points(dgat$maize_freq - dgat$teo_freq, dgat$NIROil_coeff, cex=0.5, col="red")

is.favorite_maize_enriched <- (kout_sig$maize_freq - kout_sig$teo_freq > 0 & kout_sig$NIROil_coeff > 0)
favorite_maize_enriched <- sum(is.favorite_maize_enriched)
favorite_teo_enriched <- sum(kout_sig$maize_freq - kout_sig$teo_freq < 0 & kout_sig$NIROil_coeff > 0)
unfavorite_maize_enriched <- sum(kout_sig$maize_freq - kout_sig$teo_freq > 0 & kout_sig$NIROil_coeff < 0)
is.unfavorite_teo_enriched <- (kout_sig$maize_freq - kout_sig$teo_freq < 0 & kout_sig$NIROil_coeff < 0)
unfavorite_teo_enriched <- sum(is.unfavorite_teo_enriched)
text(-0.25,0.3, labels=favorite_teo_enriched, pos=3, cex=2, col="grey")
text(-0.25,0.15, labels=paste0("(",length(dgat[,1]),")"), pos=3, cex=1.5, col="red")
text(0.25,0.3, labels=favorite_maize_enriched, pos=3, cex=2, col="darkgreen")
text(-0.25,-0.3, labels=unfavorite_teo_enriched , pos=1, cex=2, col="blue")
text(0.25,-0.3, labels=unfavorite_maize_enriched, pos=1, cex=2, col="grey")
abline(h=0, v=0, col="orange", lty=2)
dev.off()

