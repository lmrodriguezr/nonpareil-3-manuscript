#!/usr/bin/env Rscript

library(Nonpareil)
require(methods)
library(enveomics.R)

varexpl <- function(aov){
  aovss <- aov$"Sum Sq"[-1]
  cbind(aov,PctExp=c(NA,100*aovss/sum(aovss)))
}

if(FALSE){
md <- list()
for(collection in c("tara","beyster")){
  # Load metadata
  md[[collection]] <- read.table(
        paste('../OtherData/setIII-',collection,'-metadata.tsv',sep=''),
        sep='\t', header=TRUE)
  md[[collection]]$sample <- as.character(md[[collection]]$sample)
  md[[collection]]$biome <- as.factor(md[[collection]]$biome)
  rownames(md[[collection]]) <- md[[collection]]$sample

  # Estimate coverage and Nd
  md[[collection]]$nd <- NA
  md[[collection]]$mg.cov <- NA
  for(i in 1:nrow(md[[collection]])){
    f <- paste("../NPOfiles/setIII/", collection, "/",
          md[[collection]]$sample[i], ".npo", sep="")
    n <- Nonpareil.curve(f, plot=FALSE)
    md[[collection]]$nd[i] <- n$diversity
    md[[collection]]$mg.cov[i] <- n$C
  }
}

# Combine common metadata
cols <- intersect(colnames(md[[1]]), colnames(md[[2]]))
a <- rbind(md[[1]][,cols], md[[2]][,cols])

pdf("../Plots/setIII-not_included.pdf")
m <- lm(hc ~ nd, data=a)
x <- seq(11,26,.1)

col <- c("cadetblue","lightsteelblue")
plot(a$nd, a$hc, bg=col[a$biome], pch=21,
  xlab='Metagenome Nd', ylab='16S rRNA H\'', cex=(a$ssu.cov+.2)*2)

arrows(x0=a$nd, y0=a$h, y1=a$hc, col=enve.col.alpha(col[a$biome]), length=0)
p <- predict(m, a, se=TRUE, interval='prediction', level=.95)
ol <- a$hc < p$fit[,'lwr'] | a$hc > p$fit[,'upr']
text(a$nd[ol],a$hc[ol],gsub('_',' ',a$sample[ol]),col=col[a$biome[ol]],pos=1)

# Correlation model
p <- predict(m, data.frame(nd=x), se=TRUE, interval='prediction')
polygon(c(x,rev(x)), c(p$fit[,'lwr'],rev(p$fit[,'upr'])), lty=3, border='grey')

p <- predict(m, data.frame(nd=x), se=TRUE, interval='prediction',level=.9)
polygon(c(x,rev(x)), c(p$fit[,'lwr'],rev(p$fit[,'upr'])), lty=3, border='grey')

p <- predict(m, data.frame(nd=x), se=TRUE, interval='confidence')
polygon(c(x,rev(x)), c(p$fit[,'lwr'],rev(p$fit[,'upr'])), lty=2, border='grey')

abline(m, col='grey')
legend("bottomright", pt.bg=col, pch=21, legend=levels(a$biome), bty="n")
legend("topleft", paste("R:", signif(cor(a$hc, a$nd), 3)), bty="n")

dev.off()

#==> Additional data <==
a.t <- md[["tara"]]
a.b <- md[["beyster"]]
cat("Correlation\n")
cat("===========\n")
cat("## All\n")
cat("R:", cor(a$hc, a$nd), "\n")
cat("n:", nrow(a), "\n")
cat("p-value:", cor.test(a$hc, a$nd)$p.value, "\n")
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
}
cat("ANOVA Analyses\n")
cat("==============\n")
cat("## H'\n")
m <- lm(hc ~ biome + SizeFrx + LatN + Location + Wintriness + Vernality, data=a)
print(varexpl(anova(m)))
cat("## Nd\n")
m <- lm(nd ~ biome + SizeFrx + LatN + Location + Wintriness + Vernality, data=a)
print(varexpl(anova(m)))

