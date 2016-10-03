## ====================================================================
## ====================================================================
## Implements figure plotting for the following research paper:
## K. Kentzoglanakis, S. P. Brown, and R. A. Goldstein. Using the
## Price equation to analyze multi-level selection on the reproductive
## policing mechanism of bacterial plasmids
## ====================================================================
## ====================================================================


## insert all our functions
tryCatch(suppressWarnings(source('dps.r')),
         error=function(e) source('R/dps.r'))

## EPS plotting :
## setEPS(); postscript("fname.eps"); plot(); dev.off()

## ================================================================
##                 INSTRUCTIONS
## generate the pdfs and place in $PRICE_PAPER_PATH/fig.pdf/
## use gimp on each pdf and export to png (using default settings)
## place pngs in $PRICE_PAPER_PATH/figures/
## **Figure 3 ** generate fig3.pdf and fig3_inset.pdf,
##   open in gimp, paste inset into main figure and export to png
## ================================================================

## ===========================================
## === use the following files for results ===
## ===========================================

## load the results object
price.load.results <- function(fname="/data/dps/bka/results.h5") dps.load(fname)

## load the DZ object
price.load.dz <- function(fname="/data/dps/bka/dz.50k.xdr") {
  load(fname)
  dz.50k
}


## ===========================================================================
## plots the evolutionary dynamics of beta, kappa and alpha
price.plot.fig1 <- function( results, highlight=F, plot=F ) {

  if (plot)
    pdf("fig1.pdf", width=6, height=6)

  layout(matrix(1))
  ## plot plasmid replication parameters
  plot.with.range(cbind(results$dynamics$global$M$beta,
                        results$dynamics$global$M$kappa,
                        results$dynamics$global$M$alpha),
                  cbind(sqrt(results$dynamics$global$V$beta),
                        sqrt(results$dynamics$global$V$kappa),
                        sqrt(results$dynamics$global$V$alpha)),
                  col=c("blue", "red", "green"), ylim=c(0,1),
                  ylab="Plasmid Replication Parameters",
                  xlab="Evolutionary Time", xaxt='n')

  ## legend("topleft", c(expression(beta), expression(kappa), expression(alpha)),
  ##        lwd=1, col=c("blue", "red", "green"))

  if (highlight) {

    start <- 2e5
    end <- 3e5

    x <- c(start, start, end, end)
    y <- c(0, 1, 1, 0)

    polygon(x, y, border=NA, col=adjustcolor("gray", alpha.f=0.6))
    
  }

  ## mark the time ticks properly
  xticks.at <- 10000 * seq(0, 100, 25)
  axis(1, at=xticks.at)

  if (plot)
    dev.off()
  
}






## ===========================================================================
## plots the dynamics of the Price equation components for bka
price.plot.fig2 <- function(dz, steps.range=NA, ##steps.range=c(2e5, 3e5),
                            price.ylim=c(-4e-6, 4e-6), plot=F) {


  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:length(dz$fitness)
  else
    steps.range <- steps.range[1]:steps.range[2]

  first.three.quarters <- steps.range[1:(3*length(steps.range)/4)]
  last.quarter <- steps.range[(3*length(steps.range)/4):length(steps.range)]

  if (plot)
    pdf("fig2.pdf", width=10, height=5)

  layout(matrix(1:4, ncol=4), widths=c(1,rep(5,3)))
  par.bak <- par(no.readonly = TRUE)
  # plot the ylabel
  par(mar=c(6,0,4,2))
  plot(1:10, 1:10, xlab="", ylab="", type="n", axes=F)
  ##box("fig", col="red") ## for debugging
  text(5, 5.5,
       expression(paste("Price Equation Components", "  (x", 10^-6,")")),
       cex=1.5, srt=90)

  ## plot the Price equation components
  alphabet <- c("A", "B", "C")
  i <- 1

  for (name in c("beta", "kappa", "alpha")) {

    total = dz[[name]]$total[steps.range]
    inter = dz[[name]]$inter[steps.range]
    intra = dz[[name]]$intra[steps.range]
    tbias = dz[[name]]$tbias[steps.range]

    ## print out some info
    cat("==> ", name, "\n")
    cat("inter = ", mean(inter[first.three.quarters], na.rm=T),
        "|", mean(inter[last.quarter], na.rm=T), "\n")
    cat("intra = ", mean(intra[first.three.quarters], na.rm=T),
        "|", mean(intra[last.quarter], na.rm=T), "\n")

    ## plot the Price Equation components
    mplot(steps.range, cbind(total, inter, intra, tbias),
          xlab=switch(name,
              kappa=expression(paste("Evolutionary Time",
                "  (x", 10^5,")")),
              ""),
          ylab="",
          ylim=price.ylim,
          xaxt="n", yaxt="n",
          col=c("black", "blue", "red", "green"),
          cex.main=2, cex.lab=1.5)
    ## write figure label
    text(steps.range[length(steps.range)], price.ylim[2],
         labels=alphabet[i], cex=4, adj=c(1,1))
    ## draw the Y axis (major ticks)
    axis(2, at=seq(-4e-6, 4e-6, 2e-6),
         labels=switch(name, beta=seq(-4, 4, 2), NA),
         cex.axis=1.5)
    ## draw the X axis (major ticks)
    axis(1, at=seq(0, 1e6, 2.5e5),
         labels=seq(0, 10, 2.5),
         cex.axis=1.5)
    
    abline(h=0)

    i <- i + 1

  }

  par(par.bak)
  layout(matrix(1))

  if (plot)
    dev.off()

}




## ===========================================================================
## plots the dynamics of the Price equation components for bka
## for the low CN simulations
## plasmid profile params in this case are as follows:
##   phi=0.18
##   gamma=0.012
##   lambda=1
price.plot.fig2.low.cn <- function(dz, steps.range=NA,
                                   price.ylim=c(-8e-6, 8e-6), plot=F) {


  ## calculate steps.range
  if (length(steps.range) != 2)
    steps.range <- 1:length(dz$fitness)
  else
    steps.range <- steps.range[1]:steps.range[2]

  first.three.quarters <- steps.range[1:(3*length(steps.range)/4)]
  last.quarter <- steps.range[(3*length(steps.range)/4):length(steps.range)]

  if (plot)
    pdf("fig2.low.cn.pdf", width=10, height=5)

  layout(matrix(1:4, ncol=4), widths=c(1,rep(5,3)))
  par.bak <- par(no.readonly = TRUE)
  # plot the ylabel
  par(mar=c(6,0,4,2))
  plot(1:10, 1:10, xlab="", ylab="", type="n", axes=F)
  ##box("fig", col="red") ## for debugging
  text(5, 5.5,
       expression(paste("Price Equation Components", "  (x", 10^-6,")")),
       cex=1.5, srt=90)

  ## plot the Price equation components
  alphabet <- c("A", "B", "C")
  i <- 1

  for (name in c("beta", "kappa", "alpha")) {

    total = dz[[name]]$total[steps.range]
    inter = dz[[name]]$inter[steps.range]
    intra = dz[[name]]$intra[steps.range]
    tbias = dz[[name]]$tbias[steps.range]

    ## print out some info
    cat("==> ", name, "\n")
    cat("inter = ", mean(inter[first.three.quarters], na.rm=T),
        "|", mean(inter[last.quarter], na.rm=T), "\n")
    cat("intra = ", mean(intra[first.three.quarters], na.rm=T),
        "|", mean(intra[last.quarter], na.rm=T), "\n")

    ## plot the Price Equation components
    mplot(steps.range, cbind(total, inter, intra, tbias),
          xlab=switch(name,
              kappa=expression(paste("Evolutionary Time",
                "  (x", 10^5,")")),
              ""),
          ylab="",
          ylim=price.ylim,
          xaxt="n", yaxt="n",
          col=c("black", "blue", "red", "green"),
          cex.main=2, cex.lab=1.5)
    ## write figure label
    text(steps.range[length(steps.range)], price.ylim[2],
         labels=alphabet[i], cex=4, adj=c(1,1))
    ## draw the Y axis (major ticks)
    axis(2, at=seq(-8e-6, 8e-6, 2e-6),
         labels=switch(name, beta=seq(-8, 8, 2), NA),
         cex.axis=1.5)
    ## draw the X axis (major ticks)
    axis(1, at=seq(0, 1e6, 2.5e5),
         labels=seq(0, 10, 2.5),
         cex.axis=1.5)

    abline(h=0)

    i <- i + 1

  }

  par(par.bak)
  layout(matrix(1))

  if (plot)
    dev.off()

}





## ===========================================================================
## plot relatedness 
price.plot.fig3 <- function(dz, plot=F ) {

  if (plot)
    pdf("fig3.pdf", w=10, h=10)

  layout(matrix(1))

  ## plot relatedness
  mplot(Reduce(cbind,
               lapply(c("beta", "kappa", "alpha"),
                      function(z.name) dz[[z.name]]$r$cov / dz[[z.name]]$r$var),
               init=NULL),
        main="Relatedness", xlab="Evolutionary Time",
        xaxt='n')##, yaxt='n')
  ## mark the time ticks properly
  xticks.at <- 10000 * seq(0, 100, 25)
  axis(1, at=xticks.at)
  ##axis(2, at=seq(0.7,1,0.1))
      
  if (plot)
    dev.off()

}


## ===========================================================================
## plot relatedness components (cov and var)
## ==> paste into Fig3 using gimp and then convert to png
price.plot.fig3.inset <- function(dz, plot=F ) {

  if (plot)
    pdf("fig3_inset.pdf", w=4.5, h=6.25)

  par.bak <- par(no.readonly=T)

  ##layout(matrix(1:3, nrow=1))
  layout(matrix(1:2, nrow=2))

  par(mar=c(2, 4, 1, 2) + 0.1)

  ## plot covariance
  mplot(Reduce(cbind,
               lapply(c("beta", "kappa", "alpha"),
                      function(z.name) dz[[z.name]]$r$cov),
               init=NULL),
        ylab="Covariance", ##xlab="Evolutionary Time",
        log.take="y", ylim=c(-5, log10(2e-3)),
        xaxt='n')##, yaxt='n')
  ## mark the time ticks properly
  xticks.at <- 10000 * seq(0, 100, 25)
  axis(1, at=xticks.at)

  ## plot variance
  mplot(Reduce(cbind,
               lapply(c("beta", "kappa", "alpha"),
                      function(z.name) dz[[z.name]]$r$var),
               init=NULL),
        ylab="Variance", ##xlab="Evolutionary Time",
        log.take="y", ylim=c(-5, log10(2e-3)),
        xaxt='n')##, yaxt='n'
  ## mark the time ticks properly
  xticks.at <- 10000 * seq(0, 100, 25)
  axis(1, at=xticks.at)
      
  if (plot)
    dev.off()

  par(par.bak)

}



## ===========================================================================
## plot host growth and the copy number distributions over time
price.plot.fig4 <- function(results, dz, plot=F ) {

  if (plot)
    pdf("fig4.pdf", w=10, h=5)

  layout(matrix(1:2, ncol=2))

  ## plot host growth
  mplot(dz$host.growth, 
        xlab="Evolutionary Time", ylab="",
        main=expression("Host Division Rate"), xaxt='n', yaxt='n')
  ## mark the time ticks properly
  xticks.at <- 10000 * seq(0, 100, 25)
  axis(1, at=xticks.at)
  ## draw the Y axis
  axis(2, at=seq(0.05, 0.06, 0.005))
  

  ## plot copy number distributions
  x <- 0.5 + results$histograms$cn$x
  y <- t(results$histograms$cn$y)

  print(ncol(y))

  ## exclude n=0 from plot
  y <- y[2:length(x),]
  x <- x[2:length(x)]
  colors <- colorRampPalette(c("red", "blue"))(ncol(y))
  mplot(log10(x), y,
        main=expression(paste("Copy Number Frequencies (x", , 10^6, ")")),
        xlab=expression(n), 
        xaxt='n', yaxt='n', 
        col=colors)
  decorate.log("x")
  ## draw the Y axis (major ticks)
  axis(2, at=seq(0, 6e6, 2e6),
       labels=seq(0, 6, 2))
  
  abline(v=log10(7), lty="dashed")

  if (plot)
    dev.off()

}



## ===========================================================================
## plot BHS on beta with events
price.plot.fig5 <- function(dz, plot=F, steps.range=c(1,5e5)) {

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dz$beta)
  else
    steps.range <- steps.range[1]:steps.range[2]

  if (plot)
    pdf("fig5.pdf", w=10, h=4)

  layout(matrix(1:3, ncol=3))

  rng <- range(steps.range)
  xticks <- as.integer(seq(rng[1]-1, rng[2], 2.5e5))

  for (z.name in c("beta", "kappa", "alpha")) {

    nr <- dz[[z.name]]$inter.nr[steps.range]
    death <- - dz[[z.name]]$inter.death[steps.range]
    total <- nr + death
    main <- switch(z.name,
                   beta=expression(paste("BHS on ", beta, "  (x", 10^-5,")")),
                   kappa=expression(paste("BHS on ", kappa, "  (x", 10^-5,")")),
                   alpha=expression(paste("BHS on ", alpha, "  (x", 10^-5,")")))

    ## plot BHS
    mplot(steps.range, cbind(total, nr, death),
          col=c("black", "blue", "red"),
          main=main, xaxt='n', yaxt='n',
          xlab=if (z.name=="kappa") "Evolutionary Time" else "",
          ylim=c(-2e-5, 2e-5), cex.main=2, cex.lab=1.5)
    abline(h=0)
    ## mark the time ticks properly
    axis(1, at=xticks,
         labels=sapply(xticks, function(x) sprintf("%d",x)),
         cex.axis=1.5)
    ## draw the Y axis (major ticks)
    axis(2, at=seq(-2e-5, 2e-5, 2e-5),
         labels=seq(-2, 2, 2), cex.axis=1.5)
  }

  if (plot)
    dev.off()

  layout(matrix(1))

}


## ===========================================================================
## plot BHS on beta with events for low CN simulation
## plasmid profile params in this case are as follows:
##   phi=0.18
##   gamma=0.012
##   lambda=1
## This is the low-cn version of fig5
price.plot.fig6 <- function(dz, plot=F, steps.range=c(1,5e5)) {

  ## calculate steps.range
  if (length(steps.range) != 2)
    steps.range <- 1:nrow(dz$beta)
  else
    steps.range <- steps.range[1]:steps.range[2]

  if (plot)
    pdf("fig6.pdf", w=10, h=4)

  layout(matrix(1:3, ncol=3))

  rng <- range(steps.range)
  xticks <- as.integer(seq(rng[1]-1, rng[2], 2.5e5))

  for (z.name in c("beta", "kappa", "alpha")) {

    nr <- dz[[z.name]]$inter.nr[steps.range]
    death <- - dz[[z.name]]$inter.death[steps.range]
    total <- nr + death
    main <- switch(z.name,
                   beta=expression(paste("BHS on ", beta, "  (x", 10^-5,")")),
                   kappa=expression(paste("BHS on ", kappa, "  (x", 10^-5,")")),
                   alpha=expression(paste("BHS on ", alpha, "  (x", 10^-5,")")))

    ## plot BHS
    mplot(steps.range, cbind(total, nr, death),
          col=c("black", "blue", "red"),
          main=main, xaxt='n', yaxt='n',
          xlab=if (z.name=="kappa") "Evolutionary Time" else "",
          ylim=c(-4e-5, 4e-5), cex.main=2, cex.lab=1.5)
    abline(h=0)
    ## mark the time ticks properly
    axis(1, at=xticks,
         labels=sapply(xticks, function(x) sprintf("%d",x)),
         cex.axis=1.5)
    ## draw the Y axis (major ticks)
    axis(2, at=seq(-4e-5, 4e-5, 2e-5),
         labels=seq(-4, 4, 2), cex.axis=1.5)
  }

  if (plot)
    dev.off()

  layout(matrix(1))

}



## ===========================================================================
## plot plasmid parameter variance
price.plot.variance <- function(dz, plot=F ) {

  if (plot)
    pdf("variance.pdf", w=6, h=6)

  layout(matrix(1))

  y <- cbind(log10(sqrt(dz$variance$inter)),
             log10(sqrt(dz$variance$intra)))

  rng <- range(y, na.rm=T)

  mplot(y, xaxt='n', yaxt='n', ylim=c(-3, -1),
        main="Plasmid Parameter Variance", xlab="Evolutionary Time")
  decorate.log("y")
  abline(h=log10(3e-3), lty=2)

  text(10, -1, "Between Hosts", cex=2, adj=c(0, 1))
  text(10, -3, "Within Hosts", cex=2, adj=c(0, 0))

  ## mark the time ticks properly
  xticks.at <- 10000 * seq(0, 100, 25)
  axis(1, at=xticks.at)
      
  if (plot)
    dev.off()

}
