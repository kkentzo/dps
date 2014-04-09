## ====================================================================
## ====================================================================
## Implements figure plotting for the following research paper:
## K. Kentzoglanakis, S. P. Brown, and R. A. Goldstein. The evolution
## of policing in genetically mixed groups enhances productivity and
## relatedness through coercive control of neighbour reproduction.
## ====================================================================
## ====================================================================


## insert all our functions
tryCatch(suppressWarnings(source('dps.r')), error=function(e) source('R/dps.r'))

## EPS plotting :
## setEPS(); postscript("fname.eps"); plot(); dev.off()

## ========================================================
##                   INSTRUCTIONS
## generate the pdfs and place in $DPSM_PAPER_PATH/figures/
## "The invasion of cis-selfish mutants" figure is composed
## of two separate figures which should be joined using
## inkscape
## ========================================================

## === use the following files for results ===

## *** there also exist the files results.new.xdr in /data/dpsm/{htg,htg.beta}
## *** use them for the relatedness decomposition and the histograms
load.results.b <- function() dps.pp.load('/data/dpsm/htg.beta/results.xdr')
load.results.bka <- function() dps.pp.load('/data/dpsm/htg/results.xdr')
load.results.comp <- function() 
  list(ba=dps.pp.comp.load("/data/dpsm/comp/ba/results.xdr"),
       ka=dps.pp.comp.load("/data/dpsm/comp/ka/results.xdr"))

## save the results that were loaded using the above function
## in a YAML format for use with Python's plot_comp() in dps.py
## NOT IMPLEMENTED -- see plot.fig6() below
save.results.comp <- function(results.comp, fname="results.comp.yaml") {

  library(yaml)
  st <- as.yaml(results.comp)
  f <- file(fname)
  writeLines(st, f)
  close(f)
  
}


## ===========================================================================
## plots host performance (CN, growth, death)
plot.fig1 <- function(results.b, results.bka, plot=F) {

  if (plot)
    pdf("fig1.pdf", width=8, height=3)

  layout(matrix(1:3, nrow=1))

  x <- cbind(results.b$pconj, results.bka$pconj)

  ## find peak of the NO-CNC growth curve
  pconj.opt.idx <- which.max(results.b$M$custom$div.inf -
                           results.b$M$custom$death)
  ## find corresponding pconj value
  pconj.opt <- results.b$pconj[pconj.opt.idx]

  mplot(x, cbind(results.b$M$custom$cn, results.bka$M$custom$cn),
        ##cbind(results.b$M$inter$M$cn, results.bka$M$inter$M$cn),
        main="Host Copy Number", xlab=expression(p[c]), ylab="",
        log.take="x", col=c("blue", "red"))
  ##abline(v=log10(pconj.opt), h=results.b$M$custom$cn[pconj.opt.idx])
  abline(h=7, lty="dashed")
  ## legend("topleft", c("NO-CNC", "CNC"),
  ##        lwd=1, col=c("blue", "red"))

  mplot(x, cbind(results.b$M$custom$div.inf - results.b$M$custom$death,
                 results.bka$M$custom$div.inf - results.bka$M$custom$death),
        main="Host Growth", xlab=expression(p[c]), ylab="",
        log.take="x", col=c("blue", "red"))
  ##abline(v=log10(pconj.opt))

  mplot(x, cbind(results.b$M$custom$death,
                 results.bka$M$custom$death),
        main="Host Death", xlab=expression(p[c]), ylab="",
        log.take="x", col=c("blue", "red"))

  if (plot)
    dev.off()

  layout(matrix(1))
  
}





## ===========================================================================
## plots copy number distributions
plot.fig1.dists <- function(results.b.new, results.bka.new, plot=F) {

  if (plot)
    pdf("fig1dists.pdf", width=8, height=4)

  n.val <- length(results.b.new$pconj.values)

  calc.cn.hists <- function(results, take=seq(2,n.val,2)) {

    x <- 0.5 + results$histograms[[1]]$cn$x
    y <- t(Reduce(rbind,
                  lapply(results$histograms,
                         function(h) apply(h$cn$y, 2, sum)),
                  init=NULL))
    ## take the column sums for calculating probabilities
    denom <- matrix(rep(colSums(y), nrow(y)), ncol=ncol(y), byrow=T)
    list(pconj=results$pconj.values,
         x=x,
         y=(y / denom)[,take],
         m=(x %*% (y / denom))[take])
  }

  layout(matrix(1:2, nrow=1, byrow=F))

  par.bak <- par(no.readonly=T)
  
  ## plot NO-CNC histograms
  r <- calc.cn.hists(results.b.new, take=c(1,seq(2,30,2)))
  colors <- colorRampPalette(c("blue", "red"), interpolate="spline")(ncol(r$y))
  par(mar=c(4,4,3,0))
  mplot(log10(r$x), r$y,
        main="NO-CNC",
        ylab=expression(P(n)), 
        xlab=expression(n), ylim=c(0,0.12),
        xaxt='n', yaxt='n', 
        col=colors)
  decorate.log("x")
  axis(2, at=seq(0, 0.12, 0.04))
  abline(v=log10(7), lty="dashed")

  ## plot CNC histograms
  r <- calc.cn.hists(results.bka.new, take=c(1,seq(2,30,2)))
  colors <- colorRampPalette(c("blue", "red"), interpolate="spline")(ncol(r$y))
  par(mar=c(4,2,3,2))
  mplot(log10(r$x), r$y,
        main="CNC",
        xlab=expression(n), ylim=c(0,0.12),
        xaxt='n', yaxt='n', 
        col=colors, mar=c(5,3,4,2))
  decorate.log("x")
  axis(2, at=seq(0, 0.12, 0.04), labels=NA)
  abline(v=log10(7), lty="dashed")

  if (plot)
    dev.off()

  layout(matrix(1))

  par(par.bak)
  
}




## ===========================================================================
## plots plasmid performance (rep, conj, seg loss)
plot.fig2 <- function(results.b, results.bka, plot=F) {

  if (plot)
    pdf("fig2.pdf", width=8, height=8)

  layout(matrix(1:4, nrow=2, byrow=T))

  ##par.bak <- par(no.readonly=T)
  ##par(mar=0.1 + c(5, 4, 3, 2)) ## (b,l,t,r)

  ## plot rep rate
  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$rep.cn, results.bka$M$custom$rep.cn),
        main="Replication Rate", ylab="", ##xlab=expression(p[c]),,
        log.take="x")

  ## plot conj rate
  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$ht.cn, results.bka$M$custom$ht.cn),
        main="Migration Rate", ylab="", ##xlab=expression(p[c]),
        log.take="x")

  ##par(mar=0.1 + c(3, 4, 2, 2)) ## (b,l,t,r)

  ## plot seg loss
  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$loss.cn, results.bka$M$custom$loss.cn),
        main="Segregation Loss", xlab=expression(p[c]), ylab="",
        log.take="xy", col=c("blue", "red"))
  

  ## plot infection
  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$counters$inf / results.b$M$counters$n,
              results.bka$M$counters$inf / results.bka$M$counters$n),
        ##cbind(results.b$M$custom$ht.cn, results.bka$M$custom$ht.cn),
        main="Infection", ylab="", xlab=expression(p[c]),
        log.take="x")


  if (plot)
    dev.off()

  ##par(par.bak)

  layout(matrix(1))

}








## ===========================================================================
## plot the mean plasmid values and relatedness as a function of pconj
plot.fig3 <- function(results.b, results.bka, plot=F) {

  if (plot)
    pdf("fig3.pdf", width=10, height=5)

  layout(matrix(1:2, nrow=1))

  ## === plot plasmid replication parameters ===
  mplot(results.bka$pconj,
        cbind(results.bka$M$global$M$beta,
              results.bka$M$global$M$kappa,
              results.bka$M$global$M$alpha),
        main="Plasmid Parameter Values", ylim=c(0.4, 1), xlab=expression(p[c]),
        log.take="x")
  legend("left", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, box.lwd=0, col=c("blue", "red", "green"))

  ## === plot plasmid relatedness ===
  var.names <- c("beta", "kappa", "alpha")

  ## form data frames
  r.nocnc.mean <- data.frame(Reduce(cbind, results.b$M$relatedness$oo, c()))
  r.cnc.mean <- data.frame(Reduce(cbind, results.bka$M$relatedness$oo, c()))

  ## rename columns
  names(r.nocnc.mean) <- var.names
  names(r.cnc.mean) <- var.names

  mplot(results.b$pconj,
        cbind(r.nocnc.mean$beta,
              r.cnc.mean$beta,
              r.cnc.mean$kappa,
              r.cnc.mean$alpha),
        xlab=expression(p[c]), ylab="",
        ltype=c("dashed", rep("solid", 3)),
        col=c(rep("blue", 2), "red", "green"),
        log.take="x", main="Plasmid Relatedness")

  if (plot)
    dev.off()

  layout(matrix(1))

}



## this function needs the results.new.xdr files located in
## /data/dpsm/{htg,htg.beta}/
plot.fig3r <- function(results.b.new, results.bka.new, plot=F) {

  layout(matrix(1))

  if (plot)
    pdf("fig3r.pdf", w=6, h=6)

  mplot(results.b.new$pconj.values,
        cbind(
            ## plot the covariances for NO-CNC and CNC
            results.b.new$M$drelatedness$oo$cov$beta,
            results.bka.new$M$drelatedness$oo$cov$beta,
            ## plot the variances for NO-CNC and CNC
            results.b.new$M$drelatedness$oo$var$beta,
            results.bka.new$M$drelatedness$oo$var$beta
            ),
        xlab=expression(p[c]), ylab="", type="l", ltype=c(1,1,2,2),
        ylim=log10(c(3e-5, 1e-3)),
        log.take="xy", main=expression(paste("Relatedness (", beta, ")")),
        col=c("blue", "red", "blue", "red"))

  ##expression(cov(beta[ij],bar(beta)[i]))
  legend("left", c(expression(cov(beta[ij],bar(beta)[i])),
                   expression(var(beta[ij]))),
         lwd=1, col="black", lty=c(1,2))
  text(log10(results.b.new$pconj.values[15]), log10(6e-4),
       labels="CNC", cex=2, col="red")
  text(log10(results.b.new$pconj.values[15]), log10(7e-5),
       labels="NO-CNC", cex=2, col="blue")

  if (plot)
    dev.off()
  
  
}


## ===========================================================================
## plots the averages of the Price equation components as a function of pconj
## in the CNC case
plot.fig4a <- function(results.bka, plot=F) {

  if (plot)
    pdf("fig4a.pdf", w=16, h=6)

  layout(matrix(1:3, nrow=1, byrow=T))

  ## create a 1-row/3-col plot with blue:inter, red:intra and different
  ## line types for each component (n^R, n^H ...)

  for (name in c("beta", "kappa", "alpha")) {

    ylim <- switch(name,
                   beta=c(-2e-5, 2e-5),
                   kappa=c(-2e-5, 2e-5),
                   alpha=c(-2e-6, 2e-6))
    ##ylim <- NA

    main <- switch(name,
                   beta=expression(paste("Selection on ", beta, "  (x", 10^-5, ")")),
                   kappa=expression(paste("Selection on ", kappa, "  (x", 10^-5, ")")),
                   alpha=expression(paste("Selection on ", alpha, "  (x", 10^-6, ")")))

    ##ylim <- c(-2e-5, 2e-5)

    inter <- cbind(sapply(c("fitness", "nr", "ht", "death"),
                          function (fname)
                          results.bka$M$inter$C[[sprintf("%s.%s", name, fname)]],
                          USE.NAMES=F))

    intra <- cbind(sapply(c("fitness", "nr", "ht", "death"),
                          function (fname)
                          results.bka$M$intra$C[[sprintf("%s.%s", name, fname)]],
                          USE.NAMES=F))

    tbias <- results.bka$M$global$M[[sprintf("t%s", name)]]

    ## change signs in DEATH columns
    inter[,4] <- - inter[,4]
    intra[,4] <- - intra[,4]

    mplot(results.bka$pconj.values,
          cbind(inter, intra, tbias), 
          col=c(rep("blue",4), rep("red",4), "green"), ylim=ylim,
          main=main, log.take="x",
          xlab=expression(p[c]), yaxt="n",
          cex.main=2.5, cex.lab=2, cex.axis=2, ltype=1:4)

    ## draw the Y axis (major ticks)
    axis(2,
         at=switch(name,
           beta=seq(-2e-5, 2e-5, 1e-5),
           kappa=seq(-2e-5, 2e-5, 1e-5),
           alpha=seq(-2e-6, 2e-6, 1e-6)),
         labels=switch(name,
           beta=seq(-2, 2, 1),
           kappa=seq(-2, 2, 1),
           alpha=seq(-2, 2, 1)),
         cex.axis=2)
    
    abline(h=0)

    if (name == "beta") {
      legend("topleft",
             c("Selection due to f (WH)",
               expression(paste("Selection due to ", n^R, " (WH)")),
               expression(paste("Selection due to ", n^H, " (WH)")),
               expression(paste("Selection due to ", n^D, " (WH)"))),
             lwd=1, lty=1:4, bty="n", cex=2, col="red")
      legend("bottomleft",
             c("Selection due to f (BH)",
               expression(paste("Selection due to ", n^R, " (BH)")),
               expression(paste("Selection due to ", n^H, " (BH)")),
               expression(paste("Selection due to ", n^D, " (BH)"))),
             lwd=1, lty=1:4, bty="n", cex=2, col="blue")
      legend("topleft", "Transmission Bias", inset=c(0, 0.3),
             lwd=1, bty="n", cex=2, col="green")
    }
  }

  layout(matrix(1))

  if (plot)
    dev.off()

}




## ===========================================================================
## plots the averages of the Price equation components as a function of pconj
## in the NO-CNC case
plot.fig4b <- function(results.b, plot=F) {

  if (plot)
    pdf("fig4b.pdf", w=16, h=6)

  layout(matrix(1:3, nrow=1, byrow=T))

  inter <- cbind(sapply(c("fitness", "nr", "ht", "death"),
                        function (fname)
                        results.b$M$inter$C[[sprintf("beta.%s", fname)]],
                        USE.NAMES=F))
  
  intra <- cbind(sapply(c("fitness", "nr", "ht", "death"),
                        function (fname)
                        results.b$M$intra$C[[sprintf("beta.%s", fname)]],
                        USE.NAMES=F))

  tbias <- results.b$M$global$M$tbeta

  ## change signs in DEATH columns
  inter[,4] <- - inter[,4]
  intra[,4] <- - intra[,4]

  mplot(results.b$pconj.values,
        cbind(inter, intra, tbias), 
        col=c(rep("blue",4), rep("red",4), "green"),  log.take="x",
        main=expression(paste("Selection on ", beta, " (x", 10^-4,")")), 
        xlab=expression(p[c]), yaxt="n", ylim=c(-1e-4, 1e-4),
        cex.main=2.5, cex.lab=2, cex.axis=2, ltype=1:4)
  axis(2, at=seq(-1e-4, 1e-4, 1e-4), labels=c(seq(-1, 1, 1)),
       cex.axis=2)
  abline(h=0)

  legend("topleft",
         c("Selection due to f (WH)",
           expression(paste("Selection due to ", n^R, " (WH)")),
           expression(paste("Selection due to ", n^H, " (WH)")),
           expression(paste("Selection due to ", n^D, " (WH)"))),
         lwd=1, lty=1:4, bty="n", cex=2, col="red")
  legend("bottomleft",
         c("Selection due to f (BH)",
           expression(paste("Selection due to ", n^R, " (BH)")),
           expression(paste("Selection due to ", n^H, " (BH)")),
           expression(paste("Selection due to ", n^D, " (BH)"))),
         lwd=1, lty=1:4, bty="n", cex=2, col="blue")
  legend("topleft", "Transmission Bias", inset=c(0, 0.33),
         lwd=1, bty="n", cex=2, col="green")

  ## legend("topright", "Transmission Bias", ##inset=c(0, 0.35),
  ##        lwd=1, bty="n", cex=2, col="green")

  if (plot)
    dev.off()

  layout(matrix(1))
  
}




## ==================================================================================
## plot the intra-cellular covariances of (beta,alpha) and (kappa, alpha)
## as a function of pconj
plot.fig5 <- function( results, plot=F ) {

  if (plot)
    pdf("fig8.pdf", width=8, height=5)

  level = "intra"

  layout(matrix(1))

  mplot(results$pconj,
        cbind(results$M[[level]]$C$beta.alpha,
              results$M[[level]]$C$kappa.alpha),
              ##results$M[[level]]$C$alpha.fitness),
        xlab=expression(p[c]), ylab="", type="l",
        main=expression(paste("Average Within-Host Covariances",
            "  (x", 10 ^ -6, ")")),
        log.take="x", yaxt="n",
        ylim=c(-4e-6, 4e-6))
  axis(2, at=seq(-4e-6, 4e-6, 2e-6),
       labels=seq(-4, 4, 2),
       cex.axis=1)
  abline(h=0)

  ## plot error bars
  ## error.bar(log10(results$pconj),
  ##           cbind(results$M[[level]]$C$beta.alpha,
  ##                 results$M[[level]]$C$kappa.alpha),
  ##                 ##results$M[[level]]$C$alpha.fitness),
  ##           cbind(results$S[[level]]$C$beta.alpha,
  ##                 results$S[[level]]$C$kappa.alpha)  / sqrt(results$runs),
  ##           ##results$S[[level]]$C$alpha.fitness),
  ##           col=matrix(rep(c("blue", "red"),
  ##             length(results$pconj)), ncol=2, byrow=T))
  
  ## legend("topleft",
  ##        c(expression(paste("<", E[i], "[cov"[j], "(", beta, ",", alpha, ")]>")),
  ##          expression(paste("<", E[i], "[cov"[j], "(", kappa, ",", alpha, ")]>"))),
  ##        lwd=1, col=c("blue", "red"))

  if (plot)
    dev.off()
        
}





## plots the results of the comp BA and KA experiments
## produces 2 figures (fig6a, fig6b)
## open in Inkscape, paste 6b into 6a and save as pdf
plot.fig6 <- function(results.comp, plot=F) {

  ba <- results.comp$ba
  ka <- results.comp$ka
  levels <- seq(0, 0.4, by=0.01)
  nlevels <- length(levels)

  width <- 7
  height <- 6

  layout(matrix(1))
  
  ## plot BA
  if (plot)
    pdf("fig6a.pdf", w=width, h=height)
  annot <- c(0.4,0.9)
  filled.contour(ba$x.values, ba$y.values, ba$mut.wins,
                 xlab=expression(beta), ylab=expression(alpha),
                 levels=levels, nlevels=nlevels,
                 xlim=c(0.30, 0.50), ylim=c(0.8, 1),
                 col=colorpanel(length(levels), "white", "red",),##"grey10"),
                 plot.axes = {
                   axis(1); axis(2);
                   axis(3, labels=NA, tick=F); axis(4, labels=NA, tick=F)
                   points(annot[1],annot[2])
                   text(annot[1], annot[2], labels="WT", cex=1.5, pos=4)
                 } )
  if (plot)
    dev.off()
  else
    qqq <- readline("press enter to continue ")

  ## plot KA
  if (plot)
    pdf("fig6b.pdf", w=width, h=height)
  annot <- c(0.9,0.9)
  filled.contour(ka$x.values, ka$y.values, ka$mut.wins,
                 xlab=expression(kappa), ylab="",
                 levels=levels, nlevels=nlevels,
                 xlim=c(0.8, 1), ylim=c(0.8, 1),
                 col=colorpanel(length(levels), "white", "red",),##"grey10"),
                 plot.axes = {
                   axis(1); axis(2);
                   axis(3, labels=NA, tick=F); axis(4, labels=NA, tick=F)
                   points(annot[1],annot[2])
                   text(annot[1], annot[2], labels="WT", cex=1.5, pos=4)
                 } )
  if (plot)
    dev.off()

}
