#!/usr/bin/env Rscript

library(Nonpareil)
suppressMessages(library(enveomics.R))

coto.col <- c("#5BC0EB", "#FDE74C", "#9BC53D", "#E55934", 
   "#FA7921", "#EF476F", "#FFD166", "#06D6A0", "#118AB2", 
   "#073B4C", "#264653", "#2A9D8F", "#E9C46A", "#F4A261", 
   "#E76F51", "#000000")

varexpl <- function(aov){
  aovss <- aov$"Sum Sq"
  cbind(aov,PctExp=100*aovss/sum(aovss))
}

# Load data
a <- read.table('../OtherData/setI-metadata.tsv', sep='\t', header=TRUE)
a$sample <- as.character(a$sample)
a$biome <- factor(a$biome)
rownames(a) <- a$sample
a$seqtech <- factor(
      ifelse(grepl('_IonTorrent$',a$sample), "IonTorrent",
      ifelse(grepl('_454$',a$sample),"454","Illumina")))

# Estimate coverage and Nd
a$nd <- NA
a$mg.cov <- NA
for(i in 1:nrow(a)){
  f <- paste("../NPOfiles/setI/", as.character(a$biome[i]), "/",
        a$sample[i], ".npo", sep="")
  n <- Nonpareil.curve(f, plot=FALSE)
  a$nd[i] <- n$diversity
  a$mg.cov[i] <- n$C
}


m <- lm(hc ~ nd, data=a[a$h/a$hc > 0.95,])
x <- seq(11,26,.1)

pdf("../Plots/Fig_2.pdf")
col <- coto.col[ c(6,5,11,1,9,3,16) ]
plot(a$nd, a$hc, col=col[a$biome], pch=ifelse(
  grepl('_IonTorrent$',a$sample), 22, ifelse(grepl('_454$',a$sample),24,16)),
  xlab='Metagenome Nd', ylab='16S rRNA H\'', cex=(a$ssu.cov+.2)*2)
arrows(x0=a$nd, y0=a$h, y1=a$hc, col=enve.col.alpha(col[a$biome]), length=0)
p <- predict(m, a, se=TRUE, interval='prediction', level=.9)
ol <- a$hc < p$fit[,'lwr'] | a$hc > p$fit[,'upr']
text(a$nd[ol],a$hc[ol],gsub('_',' ',a$sample[ol]),col=col[a$biome[ol]],pos=1)

# Biome ranges
iqr.nd <- list()
iqr.hc <- list()
for(i in 1:length(levels(a$biome))){
  b <- levels(a$biome)[i]
  nd.q <- quantile(a$nd[a$biome==b], c(0.25,0.75), names=FALSE)
  hc.q <- quantile(a$hc[a$biome==b], c(0.25,0.75), names=FALSE)
  rect(nd.q[1], hc.q[1], nd.q[2], hc.q[2], border=enve.col.alpha(col[i]),
    col=enve.col.alpha(col[i],0.1))
  iqr.nd[[b]] <- nd.q
  iqr.hc[[b]] <- hc.q
}

# Correlation model
p <- predict(m, data.frame(nd=x), se=TRUE, interval='prediction')
polygon(c(x,rev(x)), c(p$fit[,'lwr'],rev(p$fit[,'upr'])), lty=3, border='grey')

p <- predict(m, data.frame(nd=x), se=TRUE, interval='prediction',level=.9)
polygon(c(x,rev(x)), c(p$fit[,'lwr'],rev(p$fit[,'upr'])), lty=3, border='grey')

p <- predict(m, data.frame(nd=x), se=TRUE, interval='confidence')
polygon(c(x,rev(x)), c(p$fit[,'lwr'],rev(p$fit[,'upr'])), lty=2, border='grey')

abline(m, col='grey')

legend('bottomright',
  legend=paste(gsub('-',' ',levels(a$biome)),' (n=',table(a$biome),')',sep=''),
  col=col, pch=16, bty='n', title='biome', pt.cex=1.5)
legend('topleft', legend=c('Illumina','Ion Torrent','454'), pch=c(16,22,24),
  bty='n', title='Platform')
mm <- dev.off()

pdf("../Plots/Fig_2-inset.pdf")
plot(a$ssu.cov[a$biome!="z-mock"]*100,
  (predict(m, a)-a$hc)[a$biome!="z-mock"],
  col="black", bg=col[a$biome[a$biome!="z-mock"]], pch=21,
  xlab="Turing-Good Coverage for 16S (%)",
  ylab="Residuals", cex=2)
abline(h=0)
legend("topleft", legend=cor(a$ssu.cov[a$biome!="z-mock"],
  (predict(m, a)-a$hc)[a$biome!="z-mock"]), bty="n")
mm <- dev.off()

#==> Additional data <==
a.m <-  a[a$biome=="z-mock",]
a.nm <- a[a$biome!="z-mock",]
cat("Correlation\n")
cat("===========\n")
cat("## All\n")
cat("R:", cor(a$hc, a$nd), "\n")
cat("n:", nrow(a), "\n")
cat("p-value:", cor.test(a$hc, a$nd)$p.value, "\n")
cat("## Mock\n")
cat("R:", cor(a.m$hc, a.m$nd), "\n")
cat("n:", nrow(a.m), "\n")
cat("p-value:", cor.test(a.m$hc, a.m$nd)$p.value, "\n")
cat("\n")

cat("Model Info:\n")
cat("===========\n")
print(m)
cat("\n")

cat("Difference between Nd and H\n")
cat("===========================\n")
cat("Means:", mean(a$nd)-mean(a$hc), "\n")
cat("Medians:", median(a$nd)-median(a$hc), "\n")
cat("Mean of diff:", mean(a$nd-a$hc), "\n")
cat("Median of diff:", median(a$nd-a$hc), "\n")
cat("\n")

cat("ANOVA Analyses\n")
cat("==============\n")
cat("## Simple ANOVA\n")
m <- lm(hc ~ nd, data=a)
print(varexpl(anova(m)))
cat("## Seq. Tech. ANOVA\n")
m <- lm(hc ~ nd + seqtech, data=a)
print(varexpl(anova(m)))
cat("## SSU cov ANOVA without mock\n")
m <- lm(hc ~ nd + ssu.cov, data=a.nm)
print(varexpl(anova(m)))
cat("## SSU cov ANOVA with mock\n")
m <- lm(hc ~ nd + ssu.cov, data=a)
print(varexpl(anova(m)))
cat("## MG cov ANOVA\n")
m <- lm(hc ~ nd + mg.cov, data=a.nm)
print(varexpl(anova(m)))
cat("\n")

cat("Other stats\n")
cat("===========\n")
cat("Cor. SSU-cov & MG-cov (excluding mock):",
      cor(a.nm$ssu.cov, a.nm$mg.cov), "\n")
cat("\n")

cat("IQR of Nd per biome\n")
cat("===================\n")
for(i in levels(a$biome)){
  cat(i, ":", iqr.nd[[i]], "| mid-point:", mean(iqr.nd[[i]]),"\n")
}
cat("\n")

a.a <- a[a$biome=="animal-host" & grepl("^Human_", a$sample),]
p <- quantile(a.a$nd, c(0.25,0.75))
cat("human animals: ", p, "| mid-point:", mean(p), "\n")
a.a <- a[a$biome=="animal-host" & !grepl("^Human_", a$sample),]
p <- quantile(a.a$nd, c(0.25,0.75))
cat("non-human animals: ", p, "| mid-point:", mean(p), "\n")
cat("\n")

