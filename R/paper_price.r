## insert all our functions
tryCatch(suppressWarnings(source('dps.r')),
         error=function(e) source('R/dps.r'))

## === use the following files for results ===

## || NO-CONJ results with pconj=0 || <=== THIS IS USED IN THE PAPER
## for RESULTS.BKA use results.bka <- dps.load('/data/dps/bka/results.2.16.h5')
## for DZ use load("/data/dps/bka/dz.2.16.xdr")

## EPS plotting :
## setEPS(); postscript("fname.eps"); plot(); dev.off()

## load the results to be used in the paper
price.load.bka <- function(fname="/data/dps/bka/results.15.h5") dps.load(fname)


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
                  main="", col=c("blue", "red", "green"), ylim=c(0,1),
                  xlab="Evolutionary Time", ylab="", xaxt='n')

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
price.plot.fig2 <- function(results, dz, steps.range=c(2e5, 3e5),
                            price.ylim=c(-1e-5, 1e-5), plot=F) {


  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(results$dynamics$global$M)
  else
    steps.range <- steps.range[1]:steps.range[2]

  if (plot)
    pdf("fig2.pdf", width=15, height=10)

  layout(matrix(1:6, ncol=3, byrow=F))

  alphabet <- c("A", "D", "B", "E", "C", "F")

  i <- 1

  for (name in c("beta", "kappa", "alpha")) {


    legend.loc <- switch(name,
                         beta="topright",
                         kappa="topleft",
                         alpha="topright")

    dyn <- results$dynamics$global$M[[name]][steps.range]
    dyn.sd <- sqrt(results$dynamics$global$V[[name]][steps.range])
    dyn.ylim <- range(c(dyn-dyn.sd, dyn+dyn.sd))

    ## plot the evolutionary dynamics
    plot.with.range(dyn, dyn.sd, x=steps.range, ylim=dyn.ylim,
                    main=switch(name,
                        beta=expression(paste("Evolutionary Dynamics of ", bar(beta))),
                        kappa=expression(paste("Evolutionary Dynamics of ", bar(kappa))),
                        alpha=expression(paste("Evolutionary Dynamics of ", bar(alpha)))),
                    xlab="Evolutionary Time", ylab="", col="black",
                    cex.main=2, cex.axis=1.5, cex.lab=1.5)
    
    ## write figure label
    text(steps.range[1], dyn.ylim[2], labels=alphabet[i], cex=4, adj=c(0, 1))

    i <- i+1

    total = dz[[name]]$total[steps.range]
    inter = dz[[name]]$inter[steps.range]
    intra = dz[[name]]$intra[steps.range]
    tbias = dz[[name]]$tbias[steps.range]

    ## print out some info
    cat("==> ", name, "\n")
    cat("inter = ", mean(inter, na.rm=T), "\n")
    cat("intra = ", mean(intra, na.rm=T), "\n")

    ## plot the Price Equation components
    mplot(steps.range, cbind(total, inter, intra, tbias),
          main=switch(name,
            beta=expression(paste("Price Equation Components of ",
              Delta,bar(beta), "  (x", 10^-5,")")),
            kappa=expression(paste("Price Equation Components of ",
              Delta,bar(kappa), "  (x", 10^-5,")")),
            alpha=expression(paste("Price Equation Components of ",
              Delta,bar(alpha), "  (x", 10^-5,")"))),
          ylim=price.ylim, yaxt="n",
          xlab="Evolutionary Time", ylab="", 
          col=c("black", "blue", "red", "green"),
          cex.main=2, cex.axis=1.5, cex.lab=1.5)
    ## write figure label
    text(steps.range[1], price.ylim[2], labels=alphabet[i], cex=4, adj=c(0, 1))
    ## draw the Y axis (major ticks)
    axis(2, at=seq(-1.5e-5, 1.5e-5, 0.5e-5),
         labels=seq(-1.5, 1.5, 0.5),
         cex.axis=1.5)
    ## draw the Y axis (minor ticks)
    ##axis(2, at=c(-1.5e-5, -0.5e-5, 0, 0.5e-5, 1.5e-5), labels=NA)
    
    abline(h=0)

    ## plot a frame??
    ## x <- c(rep(169466, 2), rep(177372, 2))
    ## y <- c(-1.5e-5, 1.5e-5, 1.5e-5, -1.5e-5)
    ## polygon(x, y, border=NA, col=adjustcolor("gray", alpha.f=0.6))

    i <- i + 1

  }

  if (plot)
    dev.off()

  layout(matrix(1))
  
}





## ===========================================================================
## plot the dynamics of tbiases
price.plot.fig3 <- function(dz, plot=F) {

  if (plot)
    pdf("fig3.pdf", w=6, h=6)

  layout(matrix(1))
  mplot(Reduce(cbind,
               lapply(c("beta", "kappa", "alpha"),
                      function(z.name) dz[[z.name]]$tbias),
               init=NULL), xaxt='n', yaxt='n',
        ylim=c(-3e-6, 3e-6),
        main=expression(paste("Transmission Biases (x", 10^-6,")")),
        xlab="Evolutionary Time")
  abline(h=0)
  ## mark the time ticks properly
  xticks.at <- 10000 * seq(0, 100, 25)
  axis(1, at=xticks.at)
  ## draw the Y axis (major ticks)
  axis(2, at=seq(-3e-6, 3e-6, 1.5e-6),
       labels=seq(-3, 3, 1.5))
  ##cex.axis=1.5)

  if (plot)
    dev.off()
  
  
  
}




## ===========================================================================
## plot the decline of intra beta over time
price.plot.fig4 <- function( dz, plot=F ) {

  if (plot)
    pdf("fig4.pdf", w=6, h=6)

  intra.beta <- log10(dz$beta$intra)

  rng <- range(intra.beta, na.rm=T)

  plot(log10(dz$beta$intra), t="l", xlab="Evolutionary Time", ylab="",
       main=expression(paste("Within-host Selection on ", beta)), col="blue",
       xaxt='n', yaxt='n', ylim=c(floor(rng[1]), rng[2]))
  decorate.log("y")

  ## mark the time ticks properly
  xticks.at <- 10000 * seq(0, 100, 25)
  axis(1, at=xticks.at)

  if (plot)
    dev.off()

}








## ===========================================================================
## plot relatedness 
price.plot.fig5 <- function(dz, plot=F ) {

  if (plot)
    pdf("fig5.pdf", w=6, h=6)

  layout(matrix(1))

  ## plot plasmid rep rate
  mplot(Reduce(cbind,
               lapply(c("beta", "kappa", "alpha"),
                      function(z.name) dz[[z.name]]$r),
               init=NULL),
        main="Relatedness", xlab="Evolutionary Time",
        xaxt='n', ylim=c(0.5, 1))

  ## mark the time ticks properly
  xticks.at <- 10000 * seq(0, 100, 25)
  axis(1, at=xticks.at)
      
  if (plot)
    dev.off()

}



## ===========================================================================
## plot the copy number distributions over time
price.plot.fig6 <- function(results, plot=F ) {

  if (plot)
    pdf("fig6.pdf", w=6, h=6)

  layout(matrix(1))
  
  x <- 0.5 + results$histograms$cn$x
  y <- t(results$histograms$cn$y)

  print(ncol(y))

  ## exclude n=0 from plot
  y <- y[2:length(x),]
  x <- x[2:length(x)]
  mplot(log10(x), y, main="Copy Number Distributions",
        xlab=expression(n), ylab=expression(paste("Frequency (x", 10^6, ")")),
        xaxt='n', yaxt='n',
        col=redgreen(ncol(y)))
  decorate.log("x")
  ## draw the Y axis (major ticks)
  axis(2, at=seq(0, 6e6, 1e6),
       labels=seq(0, 6, 1))
  
  abline(v=log10(7), lty="dashed")

  if (plot)
    dev.off()

}



## ===========================================================================
## plot BHS on beta with events
price.plot.fig7 <- function(dz, plot=F, steps.range=c(1,5e5)) {

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dz$beta)
  else
    steps.range <- steps.range[1]:steps.range[2]

  if (plot)
    pdf("fig6.pdf", w=10, h=6)

  layout(matrix(1:3, ncol=3))

  rng <- range(steps.range)
  xticks <- as.integer(seq(rng[1]-1, rng[2], 1e5))

  for (z.name in c("beta", "kappa", "alpha")) {

    nr <- dz[[z.name]]$inter.nr[steps.range]
    death <- - dz[[z.name]]$inter.death[steps.range]
    total <- nr + death

    main <- switch(z.name,
                   beta=expression(paste("Between-host selection on ", beta,
                     " (x", 10^-5,")")),
                   kappa=expression(paste("Between-host selection on ", kappa,
                     " (x", 10^-5,")")),
                   alpha=expression(paste("Between-host selection on ", alpha,
                     " (x", 10^-5,")")))

    ## plot BHS
    mplot(steps.range, cbind(total, nr, death),
          col=c("black", "blue", "red"),
          main=main, xaxt='n', yaxt='n',
          ylim=c(-4e-5, 4e-5), cex.main=2)
    abline(h=0)
    ## mark the time ticks properly
    axis(1, at=xticks,
         labels=sapply(xticks, function(x) sprintf("%d",x)))
    ## draw the Y axis (major ticks)
    axis(2, at=seq(-4e-5, 4e-5, 2e-5),
         labels=seq(-4, 4, 2), cex.axis=2)
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

