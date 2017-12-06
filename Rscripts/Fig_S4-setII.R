#!/usr/bin/env Rscript

library(Nonpareil)

# Load metadata
a <- read.table('../OtherData/setII-metadata.tsv', h=T, row.names=1, sep='\t')
a <- a[!is.na(a$Site), ]

# Exclude suspect contaminated datasets
excl <- read.table('../OtherData/setII-exclude.txt', h=F, sep='\t', as.is=TRUE)$V2
a <- a[!rownames(a) %in% excl, ]

# Estimate Nd
a$NpDiv <- NA
for(i in 1:nrow(a)){
  f <- paste("../NPOfiles/setII/", rownames(a)[i], ".npo", sep="")
  if(!is.na(file.size(f))){
    a$NpDiv[i] <- Nonpareil.curve(f, plot=FALSE)$diversity
  }
}

# Sort labels by NpDiv
a$Site <- factor(as.character(a$Site),
   unique(as.character(a$Site[order(a$NpDiv)])))

# Plot
pdf('../Plots/Fig_S4.pdf',9,9)
layout(1:2, heights=c(2,1.7))
par(mar=c(0,4,0.5,0.5)+0.1)
xrange <- c(13,19)# range(a$NpDiv, na.rm=TRUE)+c(-0.5,0.4)

#==> A <===
plot(1, t='n', yaxs='i', bty='l', xaxt='n', las=1, xaxs='i',
   xlab='', ylab='H\'(16S OTUs)',
   xlim=xrange, ylim=range(a$Hprime, na.rm=TRUE)*c(3.4,1.0))

legend('topright', legend='A', bty='n', text.font=2, cex=1.5)
ign <- c()
npd.m <- c()
hpr.m <- c()
mgs.c <- c()
ssu.c <- c()
n.wgs <- 0
n.ssu <- 0
col <- c("#C23499","#C2393F","#C35F3E","#C4893E","#C5C83F","#44CA5C",
      "#44C7C5","#3E57C4","#8632C5","#C492B8","#C5BC94","#95C7C6",
      rep("black",6))
for(i in 1:length(levels(a$Site))){
   color <- col[i]
   s <- levels(a$Site)[i]
   npd <- a$NpDiv[ !is.na(a$NpDiv) & a$Site == s ]
   hpr <- a$Hprime[ !is.na(a$Hprime) & a$Site == s ]
   if(length(npd)>0 & length(hpr)>0){
      npd.q <- quantile(npd, c(.05,.95), names=F)
      hpr.q <- quantile(hpr, c(.05,.95), names=F)
      points(median(npd), median(hpr), pch=16, col=color, cex=4)
      if(length(npd)>1){
	 arrows(x0=median(npd), x1=npd.q[2], y0=median(hpr), angle=90,
	    length=1/10, lwd=4, col=color)
	 arrows(x0=median(npd), x1=npd.q[1], y0=median(hpr), angle=90,
	    length=1/10, lwd=4, col=color)
      }
      if(length(hpr)>1){
	 arrows(x0=median(npd), y1=hpr.q[2], y0=median(hpr), angle=90,
	    length=1/10, lwd=4, col=color)
	 arrows(x0=median(npd), y1=hpr.q[1], y0=median(hpr), angle=90,
	    length=1/10, lwd=4, col=color)
      }
      npd.m <- c(npd.m, median(npd))
      hpr.m <- c(hpr.m, median(hpr))
      mgs.c <- c(mgs.c, length(npd))
      ssu.c <- c(ssu.c, length(hpr))
      n.wgs <- n.wgs + length(npd)
      n.ssu <- n.ssu + length(hpr)
      cat("Including: ", s, "with", length(npd), "WGS samples and",
            length(hpr), "SSU samples\n")
   }else{
      cat("Ignoring: ", s, "with", length(npd), "WGS samples and",
            length(hpr), "SSU samples\n")
      ign <- c(ign, i)
   }
}
m <- lm(hpr.m ~ npd.m)
x <- seq(13,19,.1)
y <- predict(m, data.frame(npd.m=x), interval="confidence")
polygon(c(x,rev(x)), c(y[,2], rev(y[,3])), border=NA, col="grey80")
abline(m)
legend('topleft', bty='n',
   legend=c(paste(levels(a$Site)[-ign],' (',mgs.c,', ',ssu.c,')',sep=''),
      "","","","","",""),
   col=col[-ign], ncol=2, pt.cex=2,
   pch=c(rep(16,length(levels(a$Site)[-ign])),NA,NA,NA,NA,NA,NA))
legend('bottomright', bty='n',
   legend=paste("R =", signif(cor(npd.m, hpr.m, method='pearson'),3)))
cat('n =', length(npd.m), 'sites,', n.wgs, 'WGS samples,',
      n.ssu, 'SSU samples\n')

#==> B <==
a <- a[!is.na(a$NpDiv),]
a <- a[!is.na(a$Hprime),]

par(mar=c(4,4,0.5,0.5)+0.1)
plot(1, t='n', bty='l', las=1, xaxs='i',
   xlab='Nonpareil diversity index', ylab='H\'(16S OTUs)',
   xlim=xrange, ylim=range(a$Hprime)*c(0.8,1.1))
legend('topleft', legend='B', bty='n', text.font=2, cex=1.5)
m <- lm(Hprime ~ NpDiv, data=a)
x <- seq(13,19,.1)
y <- predict(m, data.frame(NpDiv=x), interval="confidence")
polygon(c(x,rev(x)), c(y[,2], rev(y[,3])), border=NA, col="grey80")
abline(m)

points(a$NpDiv, a$Hprime, pch=16, cex=2,
   col=c('#aaaaaa55',col)[
      ifelse(is.na(a$Site), 1, as.numeric(a$Site)+1)])
legend('bottomright', bty='n',
   legend=paste("R =", signif(cor(a$NpDiv, a$Hprime, method='pearson'),3)))
cat('n =', nrow(a), '\n')

mm <- dev.off()

#==> Extra info <==

cat("Correlation\n")
cat("===========\n")
cat("R:", cor(a$NpDiv, a$Hprime), "\n")
cat("n:", nrow(a), "\n")
cat("\n")

cat("Model Info:\n")
cat("===========\n")
print(m)
cat("\n")

cat("Difference between Nd and H\n")
cat("===========================\n")
cat("Means:", mean(a$NpDiv)-mean(a$Hprime), "\n")
cat("Medians:", median(a$NpDiv)-median(a$Hprime), "\n")
cat("\n")

