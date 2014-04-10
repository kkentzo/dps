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
## *** use them for the histograms (figure 2)
load.results.b <- function(new=F) {
  if (new) dps.pp.load('/data/dpsm/htg.beta/results.new.xdr')
  else dps.pp.load('/data/dpsm/htg.beta/results.xdr')
}

load.results.bka <- function(new=F) {
  if (new) dps.pp.load('/data/dpsm/htg/results.new.xdr')
  else dps.pp.load('/data/dpsm/htg/results.xdr')
}

load.results.comp <- function() 
  list(ba=dps.pp.comp.load("/data/dpsm/comp/ba/results.xdr"),
       ka=dps.pp.comp.load("/data/dpsm/comp/ka/results.xdr"))


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

  xlab <- expression(paste("Migration rate (", p[c], ")"))

  mplot(x, cbind(results.b$M$custom$cn, results.bka$M$custom$cn),
        ##cbind(results.b$M$inter$M$cn, results.bka$M$inter$M$cn),
        main="Group Size", xlab=xlab, ylab="",
        log.take="x", col=c("blue", "red"))
  ##abline(v=log10(pconj.opt), h=results.b$M$custom$cn[pconj.opt.idx])
  abline(h=7, lty="dashed")
  ## legend("topleft", c("NO-CNC", "CNC"),
  ##        lwd=1, col=c("blue", "red"))

  mplot(x, cbind(results.b$M$custom$div.inf - results.b$M$custom$death,
                 results.bka$M$custom$div.inf - results.bka$M$custom$death),
        main="Group Performance", xlab=xlab, ylab="",
        log.take="x", col=c("blue", "red"))
  ##abline(v=log10(pconj.opt))

  mplot(x, cbind(results.b$M$custom$death,
                 results.bka$M$custom$death),
        main="Group Mortality", xlab=xlab, ylab="",
        log.take="x", col=c("blue", "red"))

  if (plot)
    dev.off()

  layout(matrix(1))
  
}





## ===========================================================================
## plots copy number distributions
## *** NEEDS THE *.new RESULTS ***
plot.fig2 <- function(results.b.new, results.bka.new, plot=F) {

  if (plot)
    pdf("fig2.pdf", width=8, height=4)

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
        main="No Policing",
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
        main="With Policing",
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
plot.fig3 <- function(results.b, results.bka, plot=F) {

  if (plot)
    pdf("fig3.pdf", width=8, height=8)

  layout(matrix(1:4, nrow=2, byrow=T))

  ##par.bak <- par(no.readonly=T)
  ##par(mar=0.1 + c(5, 4, 3, 2)) ## (b,l,t,r)

  ## plot rep rate
  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$rep.cn, results.bka$M$custom$rep.cn),
        main="Plasmid Replication Rate", ylab="", ##xlab=expression(p[c]),,
        log.take="x")

  ## plot conj rate
  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$ht.cn, results.bka$M$custom$ht.cn),
        main="Plasmid Conjugation Rate", ylab="", ##xlab=expression(p[c]),
        log.take="x")

  ##par(mar=0.1 + c(3, 4, 2, 2)) ## (b,l,t,r)

  ## plot seg loss
  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$loss.cn, results.bka$M$custom$loss.cn),
        main="Segregation Loss", ylab="",
        xlab=expression(paste("Migration rate (", p[c], ")")), 
        log.take="xy", col=c("blue", "red"))
  

  ## plot infection
  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$counters$inf / results.b$M$counters$n,
              results.bka$M$counters$inf / results.bka$M$counters$n),
        ##cbind(results.b$M$custom$ht.cn, results.bka$M$custom$ht.cn),
        main="Infection", ylab="", 
        xlab=expression(paste("Migration rate (", p[c], ")")), 
        log.take="x")


  if (plot)
    dev.off()

  ##par(par.bak)

  layout(matrix(1))

}








## ===========================================================================
## plot the mean plasmid values and relatedness as a function of pconj
plot.fig4 <- function(results.b, results.bka, plot=F) {

  if (plot)
    pdf("fig4.pdf", width=10, height=5)

  layout(matrix(1:2, nrow=1))

  ## === plot plasmid replication parameters ===
  mplot(results.bka$pconj,
        cbind(results.bka$M$global$M$beta,
              results.bka$M$global$M$kappa,
              results.bka$M$global$M$alpha),
        main="Plasmid Parameter Values", ylim=c(0.4, 1),
        xlab=expression(paste("Migration rate (", p[c], ")")),
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
        xlab=expression(paste("Migration rate (", p[c], ")")), ylab="",
        ltype=c("dashed", rep("solid", 3)),
        col=c(rep("blue", 2), "red", "green"),
        log.take="x", main="Plasmid Relatedness")

  if (plot)
    dev.off()

  layout(matrix(1))

}



## ===========================================================================
## plots the averages of the Price equation components as a function of pconj
## in the CNC case
plot.fig5a <- function(results.bka, plot=F) {

  if (plot)
    pdf("fig5a.pdf", w=16, h=6)

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
          yaxt="n", xlab=expression(p[c]),
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
plot.fig5b <- function(results.b, plot=F) {

  if (plot)
    pdf("fig5b.pdf", w=16, h=6)

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
plot.fig6 <- function( results.bka, plot=F ) {

  if (plot)
    pdf("fig6.pdf", width=8, height=5)

  level = "intra"

  layout(matrix(1))

  mplot(results.bka$pconj.values,
        cbind(results.bka$M[[level]]$C$beta.alpha,
              results.bka$M[[level]]$C$kappa.alpha),
        xlab=expression(paste("Migration rate (", p[c], ")")),
        ylab="", type="l",
        main=expression(paste("Average Within-Host Covariances",
            "  (x", 10 ^ -6, ")")),
        log.take="x", yaxt="n",
        ylim=c(-4e-6, 4e-6))
  axis(2, at=seq(-4e-6, 4e-6, 2e-6),
       labels=seq(-4, 4, 2),
       cex.axis=1)
  abline(h=0)

  if (plot)
    dev.off()
        
}




## ==================================================================================
## plots the results of the comp BA and KA experiments
## produces 2 figures (fig7a, fig7b)
## open in Inkscape, paste 7b into 7a and save as pdf
plot.fig7 <- function(results.comp, plot=F) {

  ba <- results.comp$ba
  ka <- results.comp$ka
  levels <- seq(0, 0.4, by=0.01)
  nlevels <- length(levels)

  width <- 7
  height <- 6

  layout(matrix(1))
  
  ## plot BA
  if (plot)
    pdf("fig7a.pdf", w=width, h=height)
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
    pdf("fig7b.pdf", w=width, h=height)
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
