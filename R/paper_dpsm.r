## insert all our functions
tryCatch(suppressWarnings(source('dps.r')), error=function(e) source('R/dps.r'))

## === use the following files for results ===

## || results for varying pconj ||
## for RESULTS.B use results <- dps.pp.load('/data/dpsm/htg.beta/results.xdr')
## for RESULTS.BKA use results <- dps.pp.load('/data/dpsm/htg/results.xdr')


## EPS plotting :
## setEPS(); postscript("fname.eps"); plot(); dev.off()


## ===========================================================================
## plots host performance (CN, growth, death)
plot.fig1 <- function(results.b, results.bka, plot=F) {

  if (plot)
    pdf("fig1.pdf", width=8, height=3)

  layout(matrix(1:3, nrow=1))

  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$cn, results.bka$M$custom$cn),
        main="Host Copy Number", xlab=expression(p[c]), ylab="",
        log.take="x", col=c("blue", "red"))
  legend("topleft", c("NO-CNC", "CNC"),
         lwd=1, col=c("blue", "red"))
  
  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$div.inf - results.b$M$custom$death,
              results.bka$M$custom$div.inf - results.bka$M$custom$death),
        main="Host Growth", xlab=expression(p[c]), ylab="",
        log.take="x", col=c("blue", "red"))

  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$death,
              results.bka$M$custom$death),
        main="Host Death", xlab=expression(p[c]), ylab="",
        log.take="x", col=c("blue", "red"))

  if (plot)
    dev.off()

  layout(matrix(1))
  
}






## ===========================================================================
## plots plasmid performance (rep, conj, seg loss)
plot.fig2 <- function(results.b, results.bka, plot=F) {

  if (plot)
    pdf("fig2.pdf", width=8, height=3)

  layout(matrix(1:3, nrow=1))

  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$rep.cn, results.bka$M$custom$rep.cn),
        main="Replication Rate", ylab="", xlab=expression(p[c]),
        log.take="x")
  legend("topleft", c("NO-CNC", "CNC"),
         lwd=1, col=c("blue", "red"))
  

  mplot(cbind(results.b$pconj, results.bka$pconj),
        ##cbind(results.b$M$global$M$ht, results.bka$M$global$M$ht),
        cbind(results.b$M$custom$ht.cn, results.bka$M$custom$ht.cn),
        main="Conjugation Rate", ylab="", xlab=expression(p[c]),
        log.take="x")

  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$loss.cn, results.bka$M$custom$loss.cn),
        main="Segregation Loss", xlab=expression(p[c]), ylab="",
        log.take="xy", col=c("blue", "red"))

  if (plot)
    dev.off()

  layout(matrix(1))
  
}








## ===========================================================================
## plot the mean plasmid values and relatedness as a function of pconj
plot.fig3 <- function(results.b, results.bka) {

  if (plot)
    pdf("fig3.pdf", width=8, height=6)

  layout(matrix(c(rep(1,3), 2:4), nrow=2, byrow=T))

  ## plot beta, kappa, alpha
  mplot(results.bka$pconj,
        cbind(results.bka$M$global$M$beta,
              results.bka$M$global$M$kappa,
              results.bka$M$global$M$alpha),
        main="Plasmid Parameter Values", ylim=c(0.4, 1), xlab=expression(p[c]),
        log.take="x")
  legend("left", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, box.lwd=0, col=c("blue", "red", "green"))

  ## mplot(results.b$pconj,
  ##       results.b$M$global$M$beta,
  ##       main="", xlab=expression(p[c]),
  ##       log.take="x")

  var.names <- c("beta", "kappa", "alpha")

  ## form data frames
  r.nocnc.mean <- data.frame(Reduce(cbind, results.b$M$relatedness$oo, c()))
  r.cnc.mean <- data.frame(Reduce(cbind, results.bka$M$relatedness$oo, c()))

  ## rename columns
  names(r.nocnc.mean) <- var.names
  names(r.cnc.mean) <- var.names

  for (var.name in var.names) {

    main <- switch(var.name,
                   beta=expression(paste("Relatedness (", beta, ")")),
                   kappa=expression(paste("Relatedness (", kappa, ")")),
                   alpha=expression(paste("Relatedness (", alpha, ")")))


    mplot(cbind(results.b$pconj, results.bka$pconj),
          cbind(r.nocnc.mean[[var.name]], r.cnc.mean[[var.name]]),
          xlab=expression(p[c]), ylab="", type="l",
          log.take="x", col=c("blue", "red"), main=main)

    if (var.name == "beta") {
      legend("bottomleft", c("NO-CNC", "CNC"),
             lwd=1, col=c("blue", "red"))
    }
  }

  if (plot)
    dev.off()

  layout(matrix(1))
  

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
plot.fig5 <- function( results ) {

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

