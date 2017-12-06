#!/usr/bin/env Rscript

library(Nonpareil)
suppressMessages(require(methods))
ds <- c('PFornix','Stool','Tongue','LL2011','LL2009A','LL2009B','Iowa')
kern <- c('kmer','aln95','aln99')
Nd <- list()


col <- c("#A20072","#20378C","#757575")
pdf('../Plots/Fig_S1.pdf')
layout(matrix(c(1:8,8),nrow=3, byrow=TRUE))
all <- new("Nonpareil.Set")
for(d in ds){
  Nd[[d]] <- c()
  np.set <- new("Nonpareil.Set")
  ki <- 1
  for(k in kern){
    np.out <- Nonpareil.curve(paste("../NPOfiles/div/",d,'-',k,'.npo', sep=''),
       col=col[ki], plot=FALSE)
    np.set$np.curves[[ length(np.set$np.curves)+1 ]] <- np.out 
    all$np.curves[[ length(all$np.curves)+1 ]] <- np.out 
    Nd[[d]] <- c(Nd[[d]], np.out$diversity)
    ki <- ki+1
  }
  plot(np.set, legend.opts=FALSE, plot.model=TRUE, plot.diversity=TRUE,
        arrow.head=0.04, arrow.length=0.075, main=d)
}

plot(1, t='n', xlab='Nonpareil diversity index (Nd)', xlim=range(Nd)*c(0.9,1.1),
  ylab='', ylim=c(3.5,-0.7), yaxs='i', bty='n', yaxt='n')
abline(h=1:3, col=col)
text(25, 1:3, c("Kmer","Alignment 95%","Alignment 99%"), col=col, pos=3)
for(d in ds){
  lines(Nd[[d]], 1:3, col="grey30")
  points(Nd[[d]], 1:3, col=col, pch=16, cex=3/2)
  text(Nd[[d]][1], 3/4, d, pos=4, srt=45, col="grey30")
}

mm <- dev.off()

t1 <- summary(all)
t1 <- t1[-3*(1:7),c(3,2,5)]
t1 <- cbind(t1[seq(2,14,2),], t1[seq(1,14,2),-1])[,c(1,2,4,3,5)]
t1 <- t1*rep(c(1e-9,100,100,1e-9,1e-9),each=nrow(t1))
rownames(t1) <- gsub("-.*","",rownames(t1))
colnames(t1) <- c("LR (Gbp)","C A (%)","C K (%)","LR* A (Gbp)","LR* K (Gbp)")
cat("Table 1\n")
cat("=======\n")
print(t1, digits=3)

