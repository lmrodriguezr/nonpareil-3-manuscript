#!/usr/bin/env Rscript

library(Nonpareil)
require(methods)

pdf("../Plots/Fig_S3.pdf")
layout(1:2)
col <- c(
  # 0% error
  "#EE005F",
  # 0.1% error
  "#CE005F","#CCCCCC",
  # 1% error
  "#AE005F","#888888",
  # 2% error
  "#8E005F","#333333"
  )
label.ds <- list(
      Sample1.0="Low coverage", Sample10.0="High coverage")
label.error  <- list(
      "001percent"="0.1%", onepercent="1%", twopercent="2%")
label.corr   <- list(
      error="with error correction", noerror="without error correction")
t <- list()
for(ds in c("Sample1.0", "Sample10.0")){
  k <- 1
  np.set <- new("Nonpareil.Set")
  f <- paste("../NPOfiles/err/",ds, "_zeropercent_zero.npo", sep="")
  np.set <- np.set + Nonpareil.curve(f, plot=FALSE, col=col[k],
        label="0% error")
  for(error in c("001percent","onepercent","twopercent")){
    for(corr in c("error", "noerror")){
      k <- k + 1
      f <- paste("../NPOfiles/err/", ds, "_zeropercent_",
            error, "_", corr, ".npo", sep="")
      np.set <- np.set + Nonpareil.curve(f, plot=FALSE, col=col[k],
            label=paste(label.error[[error]], label.corr[[corr]]))
    }
  }
  plot(np.set, main=label.ds[[ds]], xlim=c(1e6,1e12),
        legend.opts=list(bty="n", cex=3/4))
  t[[ length(t)+1 ]] <- summary(np.set)
}
dev.off()

cat("Supplementary Table S2\n")
cat("======================\n")
ts2 <- rbind(t(t[[1]]), t(t[[2]]))[,c("C","LRstar")]
ts2 <- cbind(ts2[c(1,2,4,6,8,9,11,13),],
      ts2[c(1,3,5,7,8,10,12,14),])[,c(1,3,2,4)]
ts2 <- ts2*rep(c(100,1e-9), each=2*nrow(ts2))
colnames(ts2) <- c("C on (%)","C off (%)","LR* on (Gbp)","LR* off (Gbp)")
rownames(ts2) <- paste(
      rep(c("0.0%","0.1%","1.0%","2.0%"),2), rep(c("LC","HC"),each=4))
print(ts2, digits=2)

