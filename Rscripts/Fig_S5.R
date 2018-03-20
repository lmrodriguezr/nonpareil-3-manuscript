#!/usr/bin/env Rscript

library(Nonpareil)
suppressMessages(library(enveomics.R))

npob <- '../NPOfiles/kmers/'
C <- list()
LRstar <- list()
Nd <- list()
K <- list()
collections <- list.files(npob)

for(collection in collections){
  f <- list.files(file.path(npob, collection), '\\.npo$')
  k <- as.numeric(gsub('.*kmer([0-9]+)\\.npo','\\1',f))
  f <- f[order(k)]
  k <- sort(k)
  K[[collection]] <- k
  C[[collection]] <- vector('numeric', length(k))
  LRstar[[collection]] <- vector('numeric', length(k))
  Nd[[collection]] <- vector('numeric', length(k))
  for(i in 1:length(k)){
    np <- Nonpareil.curve(file.path(npob, collection, f[i]), plot=FALSE)
    C[[collection]][i] <- np$C
    LRstar[[collection]][i] <- np$LRstar
    Nd[[collection]][i] <- np$diversity
  }
}

collections <- collections[ order(sapply(Nd, tail, n=1)) ]

pdf("../Plots/Fig_S5.pdf")
layout(matrix(1:4, ncol=2))
col <- rainbow(length(K), v=3/4)
names(col) <- collections

plot(1, t='n', xlim=range(K), ylim=range(C)*100, xlab='K-mer length',
      ylab='Estimated Coverage (%)', las=1)
for(i in collections) lines(K[[i]], C[[i]]*100, col=col[i], type='o', pch=16)

plot(1, t='n', xlim=range(K), ylim=range(LRstar), xlab='K-mer length',
      ylab='Projected effort (bp)', las=1, log='y')
for(i in collections) lines(K[[i]], LRstar[[i]], col=col[i], type='o', pch=16)

plot(1, t='n', xlim=range(K), ylim=range(Nd), xlab='K-mer length',
      ylab='Sequence diversity (Nd)', las=1)
for(i in collections) lines(K[[i]], Nd[[i]], col=col[i], type='o', pch=16)

plot(1, t='n', axes=FALSE, xlab='', ylab='')
legend('center', col=col, legend=collections, pch=16, bty='n')

mm <- dev.off()

