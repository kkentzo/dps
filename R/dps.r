library(tseries)
library(hdf5)
library(gplots)
library(multicore)

library(krutils)



## *** See bottom of file for command line options ***




## ============================================================================
##                     LOAD/SAVE/PROCESS HDF FILES
## ============================================================================



## load the data from file FNAME
dps.load <- function( fname ) {

  print(sprintf("Loading file %s", fname))

  results <- hdf5load(fname, load=F, tidy=T)

  results
  

}




## ============================================================================
## save the results to file FNAME
dps.save <- function( results, fname ) {

  dynamics <- results$dynamics
  settings <- results$settings
  population <- results$population
  histograms <- results$histograms

  hdf5save(fname, "dynamics", "settings", "population", "histograms")
  
}



## ============================================================================
## return a single data object containing the aggregate of all
## h5 files in the directory PATH
dps.aggregate <- function( path ) {

  fnames <- Sys.glob(file.path(path, "results.*.h5"))

  ## load the first of the files
  results <- dps.load( fnames[1] )

  for (fname in fnames[2:length(fnames)]) {

    r <- dps.load(fname)

    ## add values to dynamics
    results$dynamics$counters <- results$dynamics$counters + r$dynamics$counters

    for (level in c("global", "inter", "intra"))
      for (tp in c("M", "V", "C")) {

        results$dynamics[[level]][[tp]] <-
          results$dynamics[[level]][[tp]] +
            r$dynamics[[level]][[tp]]

      }
    
    ## add values to histograms
    results$histograms$cn$y <- results$histograms$cn$y + r$histograms$cn$y
    results$histograms$ba$z <- results$histograms$ba$z + r$histograms$ba$z

    rm(r)

  }

  ## normalize dynamics
  results$dynamics$counters <- results$dynamics$counters / length(fnames)
  
  for (level in c("global", "inter", "intra"))
    for (tp in c("M", "V", "C")) {

      results$dynamics[[level]][[tp]] <-
        results$dynamics[[level]][[tp]]  / length(fnames)
    }
  
  results
  
}











## ============================================================================
##                     PLOT RESULTS (DYNAMICS/HISTS)
## ============================================================================



## ======================================================================
## calculates and returns the correlation coefficient of the pair
## of variables specified by NAMES
## DYNAMICS should be either the global, inter or intra data frames
## STEPS.RANGE should be a 2-tuple
dps.calc.association <- function( dynamics, names, norm.by="xy", plot=F, steps.range=NA ) {

  if (length(names) == 1)
    names <- strsplit(names, "\\.")[[1]]

  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dynamics$C)
  else
    steps.range <- steps.range[1]:steps.range[2]

  ## calculate normalization factor
  norm.denom <- Reduce(`*`,
                       lapply(strsplit(norm.by, "")[[1]],
                              function (token) switch(token,
                                                      x=dynamics$V[steps.range,names[1]],
                                                      y=dynamics$V[steps.range,names[2]])),
                       1)


  assoc <- dynamics$C[steps.range,paste(names, collapse=".")] / sqrt(norm.denom)
  ##sqrt(dynamics$V[,names[1]] * dynamics$V[,names[2]])

  assoc[abs(assoc)==Inf] <- NA

  if (plot) {

    rng <- range(assoc, na.rm=T)

    x <- seq(from=rng[1], to=rng[2], by=diff(rng) / 100)
    y <- density(assoc, na.rm=T)
    ##y <- dnorm(x, mean=mean(assoc, na.rm=T), sd=sd(assoc, na.rm=T))

    h <- hist(assoc, plot=F)

    ## store settings
    par.bak <- par(no.readonly=TRUE)
    par(mfrow=c(2,1), cex.lab=1.2)

    ##title <- expression(paste("Histogram of ", rho, "(", bar(beta), ",", f, ")"))


    ## plot assoc per time frame
    breaks <- 5

    indices <- cut(steps.range, breaks=breaks, include.lowest=T, labels=F)


    colors.hist <- colorpanel(breaks, low="blue", mid="green", high="red")

    hists <- lapply(1:breaks, function (i) hist(assoc[which(indices==i)], plot=F))

    y.min <- min(sapply(hists, function (h) min(h$counts, na.rm=T)))
    y.max <- max(sapply(hists, function (h) max(h$counts, na.rm=T)))
    x.min <- min(sapply(hists, function (h) min(h$mids, na.rm=T)))
    x.max <- max(sapply(hists, function (h) max(h$mids, na.rm=T)))

    ## plot histogram
    hist(assoc, freq=F, ylim=range(h$density, y$y), xlim=c(x.min, x.max),
         xlab="Correlation Coefficient")
    ## plot density estimation
    lines(y$x, y$y, col="blue")

    means <- lapply(1:breaks, function(i) mean(assoc[which(indices==i)], na.rm=T))
    sds <- lapply(1:breaks, function(i) sd(assoc[which(indices==i)], na.rm=T))

    ## plot time-progressive histograms
    for (i in 1:breaks) {

      h <- hists[[i]]

      print(sprintf("mean=%.3f | sd=%.3f (break=%d)", means[[i]], sds[[i]], i))

      if (i==1) {

        plot(h$mids, h$counts, t="l", lwd=2, col=colors.hist[i],
             xlim=c(x.min, x.max), ylim=c(y.min, y.max))
        abline(v=0)
        
      } else {

        lines(h$mids, h$counts, lwd=2, col=colors.hist[i])
        
      }
    }

    legend("topleft", c("start", "middle", "end"),
           lwd=2, col=c("blue", "green", "red"))
    

    ## === QQ PLOT ===
    ## qqnorm(assoc, col="blue")
    ## qqline(assoc)

    ## restore settings
    par(par.bak)
    
  }

  assoc

}












## ==========================================================================
## calculates and returns the 3 components of the price equation
## for both beta and kappa
## STEPS.RANGE is a 2-tuple
dps.calc.price <- function(results, window.size=500, steps.range=NA, 
                           decompose=F, correlations=F, relatedness=T)
{

  if (length(steps.range) != 2) 
    take.steps <- 1:length(results$dynamics$global$M$fitness)
  else
    take.steps <- steps.range[1]:steps.range[2]
  
  dz <- list(fitness=ma(results$dynamics$global$M$fitness[take.steps],
                 window.size=window.size),
             window.size=window.size)

  for (z.name in c("beta", "kappa", "alpha")) {

    ## calculate and store transmission bias
    dz[[z.name]] <- data.frame(
        tbias=ma(
            results$dynamics$global$M[[sprintf("t%s", z.name)]][take.steps],
            window.size=window.size) / dz$fitness)

    ## calculate covariances
    for (level in c("global", "inter", "intra")) {

      dz[[z.name]][[level]] <-
        ma(dps.calc.association(results$dynamics[[level]],
                                sprintf("%s.fitness", z.name),
                                norm.by="",
                                steps.range=steps.range),
           window.size=window.size) / dz$fitness

      if (decompose) {

        for (name in c("nr", "ht", "death")) {

          dz[[z.name]][[sprintf("%s.%s", level, name)]] <-
            ma(dps.calc.association(results$dynamics[[level]],
                                    sprintf("%s.%s", z.name, name),
                                    norm.by="",
                                    steps.range=steps.range),
               window.size=window.size) / dz$fitness
        }
      }
    }

    ## Calculate total selection
    dz[[z.name]][["total"]] <- dz[[z.name]]$global + dz[[z.name]]$tbias

    ## Calculate relatedness
    if (relatedness) {
      dz[[z.name]]$r <- ma(results$dynamics$relatedness$wg[[z.name]][take.steps],
                           window.size=window.size)
      dz[[z.name]]$r.oo <- ma(results$dynamics$relatedness$oo[[z.name]][take.steps],
                              window.size=window.size)
    }

    
  }

  if (correlations) {

    dz$beta.alpha <- ma(dps.calc.association(results$dynamics$intra,
                                             "beta.alpha",
                                             norm.by="",
                                             steps.range=steps.range),
                        window.size=window.size)
    dz$kappa.alpha <- ma(dps.calc.association(results$dynamics$intra,
                                              "kappa.alpha",
                                              norm.by="",
                                              steps.range=steps.range),
                         window.size=window.size)
    
  }
  
  dz
  
}




## =========================================================================
## plot all the Price components (SORTED BY SELECTION LEVEL)
## DZ can be either the output of dps.calc.price() or a results objects
## use FITNESS.POSTFIX="1" to use fitness1 and "2" to use fitness2
dps.plot.price <- function( results, dz=NULL, window.size=0, steps.range=NA,
                           fitness.postfix="", xlim=NA, bw.adjust=1 ) {

  layout(rbind(rep(1,4), matrix(2:9, ncol=4)))

  ## calculate price components??
  if (is.null(dz))
    dz <- dps.calc.price(results, window.size, steps.range=NA)

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dz$beta)
  else
    steps.range <- steps.range[1]:steps.range[2]

  ## plot evolutionary dynamics of beta, kappa
  plot.with.range(cbind(results$dynamics$global$M$beta[steps.range],
                        results$dynamics$global$M$kappa[steps.range],
                        results$dynamics$global$M$alpha[steps.range]),
                  cbind(sqrt(results$dynamics$global$V$beta[steps.range]),
                        sqrt(results$dynamics$global$V$kappa[steps.range]),
                        sqrt(results$dynamics$global$V$alpha[steps.range])),
                  x=steps.range,
                  main="Plasmid Replication Parameters",
                  col=c("blue", "red", "green"),
                  xlab="Time", ylab="")
  legend("topleft", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, col=c("blue", "red", "green"))


  ## PLOT SELECTION (TOTAL, WITHIN AND BETWEEN)
  for (level in c("total", "inter", "intra", "tbias")) {

    ## plot time series
    b <- dz$beta[[level]][steps.range]
    k <- dz$kappa[[level]][steps.range]
    a <- dz$alpha[[level]][steps.range]

    print(sprintf("%s | beta=%.3e | kappa=%.3e | alpha=%.3e",
                  level, mean(b,na.rm=T),mean(k,na.rm=T), mean(a,na.rm=T)))

    mplot(steps.range, cbind(b, k, a), xlab="Time",
          main=sprintf("%s%s", level, fitness.postfix))
    abline(h=0)
    

    ## plot densities
    b <- density(b, adjust=bw.adjust, na.rm=T)
    k <- density(k, adjust=bw.adjust, na.rm=T)
    a <- density(a, adjust=bw.adjust, na.rm=T)

    mplot(cbind(b$x, k$x, a$x), cbind(b$y, k$y, a$y),
          main=sprintf("%s%s", level, fitness.postfix))
    abline(v=0)
    
  }


  print("== Discrepancies ==")
  for (name in c("beta", "kappa", "alpha")) {

    ii <- mean(dz[[name]]$inter[steps.range] + 
               dz[[name]]$intra[steps.range], na.rm=T)
    g <- mean(dz[[name]]$global[steps.range], na.rm=T)
    
    print(sprintf("%s : ABS(DIFF)=%.3e", name, abs(g-ii)))

  }

  layout(matrix(1))

}





## =========================================================================
## plot all the Price components (SORTED BY EVOLUTIONARY VARIABLES)
## DZ can be either the output of dps.calc.price() or a results objects
## use FITNESS.POSTFIX="1" to use fitness1 and "2" to use fitness2
dps.plot.price2 <- function( results, dz=NULL, window.size=0, steps.range=NA,
                            xlim=NA, bw.adjust=1 ) {

  layout(rbind(rep(1,4), matrix(2:9, ncol=4)))

  ## calculate price components??
  if (is.null(dz))
    dz <- dps.calc.price(results, window.size, steps.range=NA)

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dz$beta)
  else
    steps.range <- steps.range[1]:steps.range[2]

  ## plot evolutionary dynamics of beta, kappa, alpha
  plot.with.range(cbind(results$dynamics$global$M$beta[steps.range],
                        results$dynamics$global$M$kappa[steps.range],
                        results$dynamics$global$M$alpha[steps.range]),
                  cbind(sqrt(results$dynamics$global$V$beta[steps.range]),
                        sqrt(results$dynamics$global$V$kappa[steps.range]),
                        sqrt(results$dynamics$global$V$alpha[steps.range])),
                  x=steps.range,
                  main="Plasmid Replication Parameters",
                  col=c("blue", "red", "green"),
                  xlab="Time", ylab="")
  legend("topleft", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, col=c("blue", "red", "green"))



  ## PLOT SELECTION FOR EACH VARIABLE
  for (name in c("beta", "kappa", "alpha", "tbias")) {

    if (name == "tbias") {

      ## plot time series
      b <- dz$beta$tbias[steps.range]
      k <- dz$kappa$tbias[steps.range]
      a <- dz$alpha$tbias[steps.range]
      print(sprintf("tbias | beta=%.3e | kappa=%.3e | alpha=%.3e",
                    mean(b, na.rm=T), mean(k, na.rm=T), mean(a, na.rm=T)))
      
      mplot(steps.range, cbind(b, k, a), xlab="Time",
            main=sprintf("%s", name))
      abline(h=0)
      
      ## plot densities
      b <- density(b, adjust=bw.adjust, na.rm=T)
      k <- density(k, adjust=bw.adjust, na.rm=T)
      a <- density(a, adjust=bw.adjust, na.rm=T)

      mplot(cbind(b$x, k$x, a$x), cbind(b$y, k$y, a$y),
            main=sprintf("%s", name))
      abline(v=0)
      
      
    } else {
      ## plot time series
      total = dz[[name]]$total[steps.range]
      inter = dz[[name]]$inter[steps.range]
      intra = dz[[name]]$intra[steps.range]

      print(sprintf("%s | total=%.3e | inter=%.3e | intra=%.3e",
                    name, mean(total,na.rm=T),mean(inter,na.rm=T), mean(intra,na.rm=T)))

      mplot(steps.range, cbind(inter, intra, total), xlab="Time",
            main=sprintf("%s", name), col=c("blue", "red", "black"))
      abline(h=0)
      legend("bottomleft", c("total", "inter", "intra"),
             lwd=1, col=c("black", "blue", "red"))

      ## plot densities
      total <- density(total, adjust=bw.adjust, na.rm=T)
      inter <- density(inter, adjust=bw.adjust, na.rm=T)
      intra <- density(intra, adjust=bw.adjust, na.rm=T)

      mplot(cbind(total$x, inter$x, intra$x),
            cbind(total$y, inter$y, intra$y),
            main=sprintf("%s", name), col=c("black", "blue", "red"))
      abline(v=0)
      legend("topleft", c("total", "inter", "intra"),
             lwd=1, col=c("black", "blue", "red"))
    }
  }

  print("== Discrepancies ==")
  for (name in c("beta", "kappa", "alpha")) {

    ii <- mean(dz[[name]]$inter[steps.range] + 
               dz[[name]]$intra[steps.range], na.rm=T)
    g <- mean(dz[[name]]$global[steps.range], na.rm=T)
    
    print(sprintf("%s : ABS(DIFF)=%.3e",
                  name, abs(g-ii)))

  }

  layout(matrix(1))

}





## Use pdf("NAME.pdf", width=8, height=10) to plot to a file
dps.plot.price.name <- function(results, dz, name, steps.range=NA,
                                plot.density=F, legend.loc="topleft")
{

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(results$dynamics$global$M)
  else
    steps.range <- steps.range[1]:steps.range[2]

  layout(matrix(1:2, ncol=1))

  name.expr <- switch(name,
                      beta=expression(beta),
                      kappa=expression(kappa),
                      alpha=expression(alpha))
  dname.expr <- switch(name,
                       beta=expression(paste(Delta,beta)),
                       kappa=expression(paste(Delta,kappa)),
                       alpha=expression(paste(Delta,alpha)))
  

  ## plot the evolutionary dynamics
  plot.with.range(results$dynamics$global$M[[name]][steps.range],
                  sqrt(results$dynamics$global$V[[name]][steps.range]),
                  x=steps.range,
                  main="Evolutionary Dynamics",
                  xlab="Time", ylab="", col="black")


  total = dz[[name]]$total[steps.range]
  inter = dz[[name]]$inter[steps.range]
  intra = dz[[name]]$intra[steps.range]
  tbias = dz[[name]]$tbias[steps.range]

  ## plot the Price Equation components
  if (plot.density) {

    
    ## plot the density of Delta beta
    total <- density(dz[[name]]$total[steps.range], na.rm=T)
    inter <- density(dz[[name]]$inter[steps.range], na.rm=T)
    intra <- density(dz[[name]]$intra[steps.range], na.rm=T)
    tbias <- density(dz[[name]]$tbias[steps.range], na.rm=T)

    mplot(cbind(total$x, inter$x, intra$x, tbias$x),
          cbind(total$y, inter$y, intra$y, tbias$y),
          main="Price Equation Components", ylab="Density",
          xlab=dname.expr,
          col=c("black", "blue", "red", "green"))
    abline(v=0)

  } else {

    ## plot the dynamics
    mplot(steps.range,
          cbind(dz[[name]]$total[steps.range],
                dz[[name]]$inter[steps.range],
                dz[[name]]$intra[steps.range],
                dz[[name]]$tbias[steps.range]),
          main="Price Equation Components", xlab="Time",
          ylab="", 
          col=c("black", "blue", "red", "green"))
    abline(h=0)
    
  }
  
  legend(legend.loc, c("total", "inter", "intra", "tbias"),
         lwd=1, col=c("black", "blue", "red", "green"))

  layout(matrix(1))
  
}





dps.plot.price.prediction <- function(results, steps.range=NA) {

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(results$dynamics$global$M)
  else
    steps.range <- steps.range[1]:steps.range[2]
  

  layout(matrix(1:3, ncol=1, byrow=T))

  ## calculate price components for window size 0
  dz <- dps.calc.price(results, 0, steps.range=NA)

  ## plot dynamics
  for (name in c("beta", "kappa", "alpha")) {

    ## get actual dynamics
    adyn <- results$dynamics$global$M[[name]][steps.range]
    ## calculate one-step-ahead predicted dynamics
    pdyn <- c(adyn[1], adyn[steps.range] + dz[[name]]$total[steps.range])
    pdyn <- pdyn[1:length(pdyn)-1]

    ## plot actual and predicted
    mplot(steps.range, cbind(adyn, pdyn), main=name, xlab="Time")

    ## calculate actual dz
    dz.actual <- diff(adyn)
    d <- dz.actual[steps.range] - dz[[name]]$total[steps.range]
    print(sprintf("%s : DIFF(DZ)=%.3e", name, mean(d, na.rm=T)))
    
  }

  layout(matrix(1))
  
}






dps.plot.price.equilibrium <- function(results, dz, steps.range=NA) {

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(results$dynamics$global$M)
  else
    steps.range <- steps.range[1]:steps.range[2]

  ## calculate DZ
  ##dz <- dps.calc.price(results, 0)

  layout(matrix(1:3, nrow=1))

  for (name in c("beta", "kappa", "alpha")) {

    ## calculate densities
    total <- density(dz[[name]]$total[steps.range], na.rm=T)
    inter <- density(dz[[name]]$inter[steps.range], na.rm=T)
    intra <- density(dz[[name]]$intra[steps.range], na.rm=T)
    tbias <- density(dz[[name]]$tbias[steps.range], na.rm=T)

    if (name == "beta") {
      main <- expression(paste("Selection on ", beta))
      x <- cbind(total$x, inter$x, intra$x)
      y <- cbind(total$y, inter$y, intra$y)
      xlim <- c(-2.5e-5, 2.5e-5)
    } else if (name == "kappa") {
      main <- expression(paste("Selection on ", kappa))
      x <- cbind(total$x, inter$x, intra$x, tbias$x)
      y <- cbind(total$y, inter$y, intra$y, tbias$y)
      xlim <- c(-2e-5, 2e-5)
    } else if (name == "alpha") {
      main <- expression(paste("Selection on ", alpha))
      x <- cbind(total$x, inter$x, intra$x, tbias$x)
      y <- cbind(total$y, inter$y, intra$y, tbias$y)
      xlim <- c(-1e-5, 1e-5)
    } 

    ## plot densities
    mplot(x, y, main=main, xlim=xlim,
          col=c("black", "blue", "red", "green"), xaxt="n")
    axis(1, at=c(xlim[1], 0, xlim[2]))
    abline(v=0)
    legend("topleft", c("total", "inter", "intra", "tbias"),
           lwd=1, col=c("black", "blue", "red", "green"))
    
  }

  layout(matrix(1))
  
}





dps.plot.price.components <- function(results, dz, steps.range=NA) {

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(results$dynamics$global$M)
  else
    steps.range <- steps.range[1]:steps.range[2]

  layout(matrix(1:12, ncol=4, byrow=T))

  for (v.name in c("beta", "kappa", "alpha")) {

    for (name in c("", ".nr", ".ht", ".death")) {

      inter <- dz[[v.name]][[sprintf("inter%s", name)]][steps.range]
      intra <- dz[[v.name]][[sprintf("intra%s", name)]][steps.range]

      mplot(steps.range, cbind(inter, intra),
            main=sprintf("%s (%s)", v.name, name),
            xlab="Time", col=c("blue", "red"))
      abline(h=0)
    }
    
  }

  layout(matrix(1))
  
}








## =========================================================================
## plot the relatedness components (SORTED BY EVOLUTIONARY VARIABLES)
## DZ can be either the output of dps.calc.price() or a results objects
dps.plot.relatedness <- function( results, dz=NULL, window.size=0, 
                                 steps.range=NA, xlim=NA, bw.adjust=1 ) {

  layout(rbind(rep(1,3), matrix(2:4, ncol=3)))

  ## calculate relatedness components??
  if (is.null(dz))
    dz <- dps.calc.price(results, window.size, steps.range=NA, relatedness=T)

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dz$beta)
  else
    steps.range <- steps.range[1]:steps.range[2]

  ## plot evolutionary dynamics of beta, kappa, alpha
  plot.with.range(cbind(results$dynamics$global$M$beta[steps.range],
                        results$dynamics$global$M$kappa[steps.range],
                        results$dynamics$global$M$alpha[steps.range]),
                  cbind(sqrt(results$dynamics$global$V$beta[steps.range]),
                        sqrt(results$dynamics$global$V$kappa[steps.range]),
                        sqrt(results$dynamics$global$V$alpha[steps.range])),
                  x=steps.range,
                  main="Plasmid Replication Parameters",
                  col=c("blue", "red", "green"),
                  xlab="Time", ylab="")
  legend("topleft", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, col=c("blue", "red", "green"))



  ## PLOT RELATEDNESS FOR EACH VARIABLE
  for (name in c("beta", "kappa", "alpha")) {

    ## plot relatedness time series
    mplot(steps.range,
          cbind(dz[[name]]$r[steps.range],
                dz[[name]]$r.oo[steps.range]),
          main=sprintf("Relatedness (%s)", name),
          xlab="Time", ylab="")
    legend("bottomright", c("WG", "OO"),
           lwd=1, col=c("blue", "red"))

    print(sprintf("%s | WG=%.3e | OO=%.3e", name,
                  mean(dz[[name]]$r[steps.range] ,na.rm=T),
                  mean(dz[[name]]$r.oo[steps.range] ,na.rm=T)))

  }

  layout(matrix(1))

}



## ====================================================================================
## plot some basic features of the supplied RESULTS
dps.plot.dynamics <- function( results ) {

  ## store settings
  layout(matrix(1))

  ## plot plasmid replication parameters
  plot.with.range(cbind(results$dynamics$global$M$beta,
                        results$dynamics$global$M$kappa,
                        results$dynamics$global$M$alpha),
                  cbind(sqrt(results$dynamics$global$V$beta),
                        sqrt(results$dynamics$global$V$kappa),
                        sqrt(results$dynamics$global$V$alpha)),
                  main="Plasmid Replication Parameters",
                  col=c("blue", "red", "green"),
                  xlab="Time", ylab="")
  legend("topleft", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, col=c("blue", "red", "green"))

}



## ============================================================================
## plot an overview of RESULTS
dps.analyze <- function(results, window.size=1000) {

  ## store settings
  layout(matrix(c(1,1,2:5), ncol=2, byrow=T))


  ## plot plasmid replication parameters
  plot.with.range(cbind(results$dynamics$global$M$beta,
                        results$dynamics$global$M$kappa,
                        results$dynamics$global$M$alpha),
                  cbind(sqrt(results$dynamics$global$V$beta),
                        sqrt(results$dynamics$global$V$kappa),
                        sqrt(results$dynamics$global$V$alpha)),
                  main="Plasmid Replication Parameters",
                  col=c("blue", "red", "green"),
                  xlab="Time", ylab="")
  legend("topleft", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, col=c("blue", "red", "green"))


  ## plot division rate
  div.rate <- ma(results$dynamics$counters$div.inf / results$dynamics$counters$n,
                 window.size=window.size)
  plot(div.rate, col="blue", t="l", xlab="Time", ylab="",
       main="Average Host Division Rate")


  ## plot ave copy number
  cn <- ma(results$dynamics$counters$cn / results$dynamics$counters$n, window.size=window.size)
  plot(cn, col="blue", t="l", xlab="Time", ylab="",
       main="Average Copy Number")



  ## plot replication rate per plasmid
  rep.rate <- ma(results$dynamics$counters$rep / results$dynamics$counters$cn, window.size=window.size)
  plot(rep.rate, col="blue", t="l", xlab="Time", ylab="",
       main="Average Plasmid Replication Rate")


  ## plot replication rate per plasmid
  ht.rate <- ma(results$dynamics$counters$ht / results$dynamics$counters$cn, window.size=window.size)
  plot(ht.rate, col="blue", t="l", xlab="Time", ylab="",
       main="Average HT Rate")

  layout(matrix(1))
  
  
}








dps.plot.histogram <- function( results, name="age" ) {

  x <- results$histograms[[name]]$x
  y <- colSums(results$histograms[[name]]$y)

  x <- x[y!=0]
  y <- y[y!=0]

  plot(x, log10(y), t="b", xlab=name, ylab="Frequency", yaxt="n")
  decorate.log("y")
  
  
}






## ============================================================================
## CMD LINE FUNCTION : calls dps.analyze() and dps.plot.price() for all files in H5.PATH and
## saves the generated png files in PNG.PATH
## **multicore function : specify CORES **
dps.plot.all.parallel <- function( h5.path, png.path, cores ) {


  ## get the file names of all the results 
  fnames <- Sys.glob(file.path(h5.path, "results.[0-9]*.h5"))

  ## are there any results in path??
  if (length(fnames) == 0) 
    stop(paste("No results found in supplied path:", h5.path))

  ## if the png dir does not exist - create it
  if (! file.exists(png.path))
    dir.create(png.path, recursive=T)


  plot.results <- function( fnames ) {

    for (fname in fnames) {

      ## split filename tokens
      tokens <- strsplit(basename(fname), "\\.")[[1]]
      ## form for generating figure file names
      r.id <- paste(tokens[2], tokens[3], sep=".")

      ## load results
      r <- dps.load(fname)

      ## analyze results
      png(file.path(png.path, paste("results", r.id, "png", sep=".")),
          width=640, height=640, bg="white")
      dps.plot.dynamics(r)
      ##dps.analyze(r, window.size=500)
      dev.off()

      ## calculate DZ
      ##dz <- dps.calc.price(r, 500)

      ## save DZ
      ##save(dz, file=file.path(png.path, paste("dz", r.id, "xdr", sep=".")), compress=T)

      ## plot price dynamics
      ## png(file.path(png.path, paste("price", r.id, "png", sep=".")),
      ##     width=640, height=640, bg="white")
      ## dps.plot.price(r, dz, 0)
      ## dev.off()

      ## rm(dz)
      rm(r)

      print(sprintf("DONE : %s", r.id))
      
    }

    fnames
    
  }

  ## =====================================================================================
  ## run the plotting in parallel
  pvec(fnames, plot.results, mc.cores=cores)

}







## ============================================================================
## CMD LINE FUNCTION : calls dps.analyze() and dps.plot.price() for
## the results (h5) file specified as FNAME
## saves the generated png (and xdr) files in PNG.PATH
dps.plot.all.single <- function( fname, png.path, cores ) {

  ## if the png dir does not exist - create it
  if (! file.exists(png.path))
    dir.create(png.path, recursive=T)

  ## split filename tokens
  tokens <- strsplit(basename(fname), "\\.")[[1]]
  ## form for generating figure file names
  r.id <- paste(tokens[2], tokens[3], sep=".")

  ## load results
  r <- dps.load(fname)

  ## analyze results
  ## png(file.path(png.path, paste("results", r.id, "png", sep=".")),
  ##     width=640, height=640, bg="white")
  ## dps.analyze(r, window.size=500)
  ## dev.off()

  ## calculate DZ
  print("Calculating DZ")
  dz <- dps.calc.price(r, 1000)

  ## save DZ
  print(file.path(png.path, paste("dz", r.id, "xdr", sep=".")))
  save(dz, file=file.path(png.path, paste("dz", r.id, "xdr", sep=".")), compress=T)

  ## plot price dynamics for beta, kappa
  print("Plotting Price dynamics")
  png(file.path(png.path, paste("price", r.id, "png", sep=".")),
      width=640, height=640, bg="white")
  dps.plot.price(r, dz, 0)
  dev.off()

  print(sprintf("DONE : %s", r.id))

}








## ============================================================================
## CMD LINE FUNCTION : called from the command line
## if CORES > 1, it calls dps.plot.all.parallel()
## otherwise, it calls dps.plot.all.single()
dps.plot.all <- function( h5.path, png.path, cores ) {

  if (cores > 1)
    dps.plot.all.parallel(h5.path, png.path, cores)
  else
    dps.plot.all.single(h5.path, png.path, 1)
  
}




## ============================================================================
## plots the joint histogram of the requested pair
dps.plot.joint.dist <- function( results, pair="bk", fparts=NA, lv=15 ) {


  get.expression <- function( axis ) {

    switch(substr(pair, axis, axis), b=expression(beta), k=expression(kappa),
           a=expression(alpha), NA)
    
  }

  ## is results a path to results??
  if (is.character(results)) 
    ## load results from that path
    results <- dsps.load(results)

  x <- results$histograms[[pair]]$x
  y <- results$histograms[[pair]]$y

  if (all(is.na(fparts)))
    fparts <- 1:dim(results$histograms[[pair]]$z)[1]
  
  steps.per.fpart <- results$settings$steps %/% dim(results$histograms[[pair]]$z)[1]
  n.rows <- floor(sqrt(length(fparts)))
  n.cols <- length(fparts) %/% n.rows + ifelse(length(fparts) %% n.rows, 1, 0)

  layout(matrix(1:(n.rows*n.cols), nrow=n.rows, ncol=n.cols, byrow=T))

  levels <- pretty(results$histograms[[pair]]$z[fparts,,], lv)

  for (fpart in fparts) {
    
    z <- results$histograms[[pair]]$z[fpart,,] / sum(results$histograms[[pair]]$z[fpart,,])

    image(x, y, z, col=colorpanel(length(levels), low="white", high="red"),
          ##col=colorpanel(length(levels), low="white", mid="grey", high="black"),
          xlab=get.expression(1), ylab=get.expression(2), 
          main=sprintf("Steps %d-%d", (fpart-1)*steps.per.fpart, fpart*steps.per.fpart))
    
    abline(v=x, lty="dotted", lwd=0.5)
    abline(h=y, lty="dotted", lwd=0.5)

  }

  layout(matrix(1))

}







### ===================================================================================
### ===================================================================================
###                      CLUSTER EXPERIMENTS
### ===================================================================================
### ===================================================================================





## load the XDR file with the post-processed results
dps.pp.load <- function( fname ) {

  load(fname)
  results
  
}








## =====================================================================================
## CMD LINE FUNCTION :: post-process the results of the HT.GLOBAL experiment concurrently
## pass NA to cores and it will just print the dimensionality (pconj,runs) and then exit
dps.pp.experiment.parallel <- function( path, cores ) {

  ## extract the first number from the string FNAME
  ## this number is the id for the index of pconj
  ## return this id
  extract.id <- function( fname, extract.digit=c("pconj", "run") ) {

    if (length(extract.digit) > 1)
      stop("please provide a single string to extract.digit")

    tokens <- strsplit(fname, "\\.")[[1]]

    pos <- switch(extract.digit, pconj=2, run=1, NA)

    as.numeric(tokens[length(tokens) - pos])
    
  }

  ## get the file names of all the results 
  fnames <- Sys.glob(file.path(path, "results.[0-9]*.h5"))

  ## are there any results in path??
  if (length(fnames) == 0) 
    stop(paste("No results found in supplied path:", path))

  ## count number of gamma, cn, sigma and runs values
  len.pconj <- 0
  runs <- 0
  for (fname in fnames) {
    
    ## extract values
    pconj.value <- extract.id(fname, "pconj")
    run.value <- extract.id(fname, "run")

    ##print(sprintf("%s (%d,%d)", fname, pconj.value, run.value))

    if (run.value == 0)
      len.pconj <- len.pconj + 1
    if (pconj.value == 0)
      runs <- runs + 1
  }

  print(sprintf("Dimensionality=(pconj=%d,runs=%d)", len.pconj, runs))

  ## having printed the dim info, exit or continue??
  if (is.na(cores))
    return()
  else
    print(sprintf("Number of cores to use : %d", cores))



  ## =====================================================================================
  ## process the results for the specified index of pconj (for all runs)
  ## return a RESULTS list
  process.pconj <- function( i.pconj ) {

    pconj.values <- NA

    results <- list(counters=list(), ## automatic
                    relatedness=list(), ## automatic
                    global=list(), ## automatic
                    inter=list(), ## automatic
                    intra=list(), ## automatic
                    custom=list(
                        div.inf=array(NA, runs),
                        div.all=array(NA, runs),
                        death=array(NA, runs),
                        extinction=array(0, runs),
                        cn=array(NA, runs),
                        rep.cn=array(NA, runs),
                        rep.n=array(NA, runs),
                        ht.cn=array(NA, runs),
                        ht.n=array(NA, runs),
                        loss.cn=array(NA, runs),
                        loss.n=array(NA, runs),
                        inf=array(NA, runs)))

    
    for (i.run in 1:runs) {
      ## form file name
      fname <- sprintf("results.%d.%d.h5", i.pconj-1, i.run-1)

      ## load results
      r <- dps.load(file.path(path, fname))

      ## get a sequence of steps indexing the second half of the simulation
      steps <- nrow(r$dynamics$global$M)
      seq.steps <- (ceiling(steps/5)):steps

      ## add value of gamma to results (as an atomic value)
      pconj.values <- r$settings$pconj

      ## do we have extinction?
      if (steps < r$settings$steps && r$dynamics$counters$ptypes[steps] == 0) 
        results$custom$extinction[i.run] <- 1

      ## process CUSTOM elements
      results$custom$div.inf[i.run] <- mean(r$dynamics$counters$div.inf[seq.steps] /
                                            r$dynamics$counters$n[seq.steps],
                                            na.rm=T)
      results$custom$div.all[i.run] <- mean(r$dynamics$counters$div.all[seq.steps] /
                                            r$dynamics$counters$n[seq.steps],
                                            na.rm=T)
      results$custom$death[i.run] <- mean(r$dynamics$counters$death[seq.steps]
                                          / r$dynamics$counters$n[seq.steps], na.rm=T)
      results$custom$cn[i.run] <- mean(r$dynamics$counters$cn[seq.steps]
                                       / r$dynamics$counters$n[seq.steps], na.rm=T)
      results$custom$rep.cn[i.run] <- mean(r$dynamics$counters$rep[seq.steps]
                                           / r$dynamics$counters$cn[seq.steps], na.rm=T)
      results$custom$rep.n[i.run] <- mean(r$dynamics$counters$rep[seq.steps]
                                          / r$dynamics$counters$n[seq.steps], na.rm=T)
      results$custom$ht.cn[i.run] <- mean(r$dynamics$counters$ht[seq.steps]
                                          / r$dynamics$counters$cn[seq.steps], na.rm=T)
      results$custom$ht.n[i.run] <- mean(r$dynamics$counters$ht[seq.steps]
                                         / r$dynamics$counters$n[seq.steps], na.rm=T)
      results$custom$loss.cn[i.run] <- mean(r$dynamics$counters$loss[seq.steps]
                                            / r$dynamics$counters$cn[seq.steps], na.rm=T)
      results$custom$loss.n[i.run] <- mean(r$dynamics$counters$loss[seq.steps]
                                           / r$dynamics$counters$n[seq.steps], na.rm=T)
      results$custom$inf[i.run] <- mean(r$dynamics$counters$inf[seq.steps]
                                        / r$dynamics$counters$n[seq.steps], na.rm=T)

      ## process COUNTERS
      for (name in names(r$dynamics$counters)) {
        ## initialize array (for i.run==1)??
        if (is.null(results$counters[[name]])) 
          results$counters[[name]] <- array(NA, runs)

        ## record mean
        results$counters[[name]][i.run] <- mean(r$dynamics$counters[[name]][seq.steps],
                                                na.rm=T)
        
      }

      ## process RELATEDNESS
      for (r.type in names(r$dynamics$relatedness)) 
        for (name in names(r$dynamics$relatedness[[r.type]])) {
          if (is.null(results$relatedness[[r.type]][[name]])) 
            results$relatedness[[r.type]][[name]] <- array(NA, runs)
          
          results$relatedness[[r.type]][[name]][i.run] <-
            mean(r$dynamics$relatedness[[r.type]][[name]][seq.steps], na.rm=T)
          
        }
      

      ## process the rest
      for (level in c("global", "inter", "intra"))
        for (stat in names(r$dynamics[[level]]))
          for (name in names(r$dynamics[[level]][[stat]])) {
            ## initialize array??
            if (is.null(results[[level]][[stat]][[name]])) {
              results[[level]][[stat]][[name]] <- array(NA, runs)
              ## intialize array for correlation coefficients
              if (stat == "C")
                results[[level]][["CC"]][[name]] <- array(NA, runs)
            }

            ## record average statistic
            results[[level]][[stat]][[name]][i.run] <-
              mean(r$dynamics[[level]][[stat]][[name]][seq.steps], na.rm=T)
            ## calculate and record correlation coefficient in CC
            if (stat == "C")
              results[[level]][["CC"]][[name]][i.run] <-
                mean(dps.calc.association(r$dynamics[[level]],
                                          name,
                                          plot=F)[seq.steps],
                     na.rm=T)
          }

      ## remove results from memory
      rm(r)
      
    }

    ## form aggregated results
    results.agg <- list(M=list(counters=list(),relatedness=list(), custom=list()),
                        S=list(counters=list(), relatedness=list(), custom=list()),
                        pconj.values=pconj.values)

    ## aggregate the results into results.agg
    for (level in names(results)) {

      if (level == "counters" || level == "custom") {

        for (name in names(results[[level]])) {
          results.agg$M[[level]][[name]] <- mean(results[[level]][[name]], na.rm=T)
          results.agg$S[[level]][[name]] <- sd(results[[level]][[name]], na.rm=T)
        }

        
      } else if (level == "relatedness") {

        for  (r.type in names(results[[level]])) {

          if (is.null(results.agg$M[[level]][[r.type]])) {
            results.agg$M[[level]][[r.type]] <- list()
            results.agg$S[[level]][[r.type]] <- list()
          }
          
          for (name in names(results[[level]][[r.type]])) {

            results.agg$M[[level]][[r.type]][[name]] <-
              mean(results[[level]][[r.type]][[name]], na.rm=T)
            results.agg$S[[level]][[r.type]][[name]] <-
              sd(results[[level]][[r.type]][[name]], na.rm=T)
          }
        }
        
      } else{

        for (stat in names(results[[level]])) {
          if (is.null(results.agg$M[[level]][[stat]])) {
            results.agg$M[[level]][[stat]] <- list()
            results.agg$S[[level]][[stat]] <- list()
          }

          for (name in names(results[[level]][[stat]])) {
            results.agg$M[[level]][[stat]][[name]] <-
              mean(results[[level]][[stat]][[name]], na.rm=T)
            results.agg$S[[level]][[stat]][[name]] <-
              sd(results[[level]][[stat]][[name]], na.rm=T)
          }
        }
      }
    }

    results.agg
    
  }  
  



  ## =====================================================================================
  ## AD-HOC TESTING (COMMENT OUT WHEN DONE)
  ## len.pconj <- 1
  ## results <- process.pconj( 1 )
  ## return(results)
  ## =====================================================================================

  
  
  ## =====================================================================================
  ## split the i_pconj.* processing according to the cores we have
  indices <- cut(1:len.pconj, breaks=seq(0, len.pconj+cores, by=cores),
                 include.lowest=T, labels=F)

  indices.rng <- range(indices)

  results <- list()

  for (i in indices.rng[1]:indices.rng[2]) {

    ## assemble the values of pconj indices for this batch
    i.pconj.values <- (1:len.pconj)[which(indices == i)]
    ## run jobs on this batch
    jobs <- lapply(i.pconj.values, function(x) parallel(process.pconj(x), name=x))
    results <- c(results, collect(jobs, wait=TRUE))

  }


  ## do we have only one element in the list???
  if (len.pconj > 1) {

    ## === AGGREGATE RESULTS ===

    for (i.pconj in 2:len.pconj) {
      ## print value of i.pconj
      print(sprintf("i=%d", i.pconj))
      ## first, the PCONJ.VALUES
      results[["1"]]$pconj.values <- c(results[["1"]]$pconj.values,
                                       results[[as.character(i.pconj)]]$pconj.values)
      ## second, the LEVELS (counters, global, inter, intra, relatedness)
      for (level in names(results[["1"]]$M)) {
        if (level == "counters" || level == "custom") {
          for (name in names(results[["1"]]$M[[level]])) {
            results[["1"]]$M[[level]][[name]] <-
              c(results[["1"]]$M[[level]][[name]],
                results[[as.character(i.pconj)]]$M[[level]][[name]])
            results[["1"]]$S[[level]][[name]] <-
              c(results[["1"]]$S[[level]][[name]],
                results[[as.character(i.pconj)]]$S[[level]][[name]])
          }
        } else if (level == "relatedness") {

          for (r.type in names(results[["1"]]$M[[level]]))
            for (name in names(results[["1"]]$M[[level]][[r.type]])) {

              results[["1"]]$M[[level]][[r.type]][[name]] <-
                c(results[["1"]]$M[[level]][[r.type]][[name]],
                  results[[as.character(i.pconj)]]$M[[level]][[r.type]][[name]])
              results[["1"]]$S[[level]][[r.type]][[name]] <-
                c(results[["1"]]$S[[level]][[r.type]][[name]],
                  results[[as.character(i.pconj)]]$S[[level]][[r.type]][[name]])
            }
          
        } else {

          for (stat in names(results[["1"]]$M[[level]]))
            for (name in names(results[["1"]]$M[[level]][[stat]])) {
              results[["1"]]$M[[level]][[stat]][[name]] <-
                c(results[["1"]]$M[[level]][[stat]][[name]],
                  results[[as.character(i.pconj)]]$M[[level]][[stat]][[name]])
              results[["1"]]$S[[level]][[stat]][[name]] <-
                c(results[["1"]]$S[[level]][[stat]][[name]],
                  results[[as.character(i.pconj)]]$S[[level]][[stat]][[name]])
              
            }
        }
      }
    }
  }
  
  ## form results
  results <- results[["1"]]

  ## store number of runs
  results$runs <- runs

  ## save as HDR 
  save(results, file=file.path(path, "results.xdr"), compress=T)

}









## ======================================================================================
## plot a variable from the PP results, specified by TYPE, LEVEL and NAME
## use repeatedly on the same figure using ADD=TRUE
## if PLOT.ERRORS, then error bars with standard error of the mean will be plotted 
dps.pp.plot <- function(results, type="M", level="inter", name="bk",
                        main=NA, ylab="", ylim=NA, col="blue",
                        plot.errors=F, add=F, take.log="")
{

  if (is.na(as.character(main)))
    main <- sprintf("%s (%s)", name, level)

  if (level == "custom" || level == "counters") {
    y <- results$M[[level]][[name]]
    err <- results$S[[level]][[name]]
  } else {
    y <- results$M[[level]][[type]][[name]]
    err <- results$S[[level]][[type]][[name]]
  }

  ## calculate standard error
  err <- err / sqrt(results$runs)

  if (length(ylim) != 2)
    ylim <- range(c(y-err, y+err), na.rm=T)

  if (add)
    lines(log10(results$pconj.values), y, col=col)
  else
    plot(log10(results$pconj.values), y, t="l",
         ylim=ylim,
         xlab=expression(p[c]), xaxt='n', col=col, ylab="",
         main=main)

  ## plot error bars
  if (plot.errors)
    error.bar(log10(results$pconj.values), y, err,
              length=0.05, col=col)

  ## draw logarithmic axis??
  if (! add)
    decorate.log("x")
  
}








dps.pp.plot.price <- function(results) {

  layout(matrix(1:6, ncol=3))

  x <- results$pconj.values

  for (level in c("global", "inter", "intra")) {

    ##for (fitness.postfix in c("", "1", "2")) {

    fitness.postfix <- ""

    y <- cbind(results$M[[level]]$C[[sprintf("beta.fitness%s", fitness.postfix)]],
               results$M[[level]]$C[[sprintf("kappa.fitness%s", fitness.postfix)]],
               results$M[[level]]$C[[sprintf("alpha.fitness%s", fitness.postfix)]])
    
    mplot(x, y, main=sprintf("%s selection %s", level, fitness.postfix), log.take="x")
    decorate.log("x")
    abline(h=0)
    ##}
  }

  ## plot total selection
  mplot(x, cbind(results$M$inter$C$beta.fitness + results$M$intra$C$beta.fitness
                 ##results$M$global$C$beta.fitness
                 + results$M$global$M$tbeta,
                 results$M$inter$C$kappa.fitness + results$M$intra$C$kappa.fitness
                 ##results$M$global$C$kappa.fitness
                 + results$M$global$M$tkappa,
                 results$M$inter$C$alpha.fitness + results$M$intra$C$alpha.fitness
                 ##results$M$global$C$alpha.fitness
                 + results$M$global$M$talpha),
        main="Total Selection", log.take="x")
  decorate.log("x")
  abline(h=0)


  ## plot transmission bias
  mplot(x, cbind(results$M$global$M$tbeta, results$M$global$M$tkappa,
                 results$M$global$M$talpha),
        main="Transmission Bias", log.take="x")
  decorate.log("x")
  abline(h=0)
  
  ## plot discrepancy
  mplot(results$pconj.values,
        cbind(results$M$global$C$beta.fitness -
              (results$M$inter$C$beta.fitness + results$M$intra$C$beta.fitness),
              results$M$global$C$kappa.fitness -
              (results$M$inter$C$kappa.fitness + results$M$intra$C$kappa.fitness),
              results$M$global$C$alpha.fitness -
              (results$M$inter$C$alpha.fitness + results$M$intra$C$alpha.fitness)),
        main="Discrepancy", log.take="x")
  decorate.log("x")
  abline(h=0)
  
  layout(matrix(1))
  
}





dps.pp.plot.price.components <- function(results, name) {

  layout(matrix(1:4, nrow=1, byrow=T))

  for (fname in c("fitness", "nr", "ht", "death")) {

    inter <- results$M$inter$C[[sprintf("%s.%s", name, fname)]]
    intra <- results$M$intra$C[[sprintf("%s.%s", name, fname)]]

    if (fname == "nr") {
      main <- expression(paste("Selection due to ", n^R))
    } else if (fname == "ht") {
      main <- expression(paste("Selection due to ", n^H))
    } else if (fname == "death") {
      main <- expression(paste("Selection due to ", n^D))
      inter = - inter
      intra = - intra
    } else if (fname == "fitness") {
      main <- "Total Selection"
    }

    print(sprintf("INTER | %s = %.3e", fname, inter[length(inter)]))
    print(sprintf("INTRA | %s = %.3e", fname, intra[length(intra)]))

    mplot(log10(results$pconj.values),
          cbind(inter, intra), 
          xaxt='n', col=c("blue", "red"),
          main=main,
          xlab=expression(p[c]))
    abline(h=0)
    decorate.log("x")

    if (fname == "fitness")
      legend("topleft", c("inter", "intra"),
             lwd=1, col=c("blue", "red"))

  }

  layout(matrix(1))
  
}





## ======================================================================================
## ======================================================================================
##                             CREATE CLUSTER SCRIPTS
## ======================================================================================
## ======================================================================================


## Generates a script for the `sbatch` cluster command
## the LABEL should also be a directory in /home/kkentzo/
## using an sprintf command to substitute the value which will be determined
## according to the contents of LOAD.POP.NAMES
## `...` contains arguments to dsps:
##        use . instead of _ (e.g. ht.model) and it will converted automatically to _

## this is a generic function --
## see dps.create.experiment.script() below for more specific experiments
dps.create.cluster.script <- function(label="htg", fname=NA, cores=100, runs=10,
                                      pconj.values=seq(0.02, 0.15, length.out=50),
                                      ...) {

  ## adjust number of cores??
  total.no.of.jobs <- length(pconj.values) * runs
  if (total.no.of.jobs < cores)
    cores <- total.no.of.jobs
  

  ## start with header
  contents <- paste("#!/bin/sh\n",
                    sprintf("#SBATCH --job-name=%s", toupper(label)),
                    sprintf("#SBATCH --output=%s.out", label),
                    sprintf("#SBATCH --error=%s.err", label),
                    sprintf("#SBATCH --ntasks=%d", cores),
                    "\ndate\n\n", sep="\n")

  dps.path <- "LD_LIBRARY_PATH=/home/kkentzo/local/lib /home/kkentzo/Projects/programs/dps"


  ## construct extra arguments string
  args.list <- list(...)

  ## convert steps string from %e to %d format
  if (length(args.list) > 0 && ! is.na(args.list$steps)) {
    args.list$steps <- sprintf("%d", args.list$steps)
  }

  ## assemble all arguments in one string using the [--arg.name arg.value] format
  args.string <- paste(lapply(names(args.list),
                              function (arg.name) {
                                sprintf("--%s %s", paste(strsplit(arg.name, "\\.")[[1]],
                                                         collapse="_"),
                                        args.list[[arg.name]])
                              }),
                       collapse=" ")

  ## create commands
  job.id <- 0
  wait.inserted <- F
  
  for (i.pconj in 0:(length(pconj.values)-1))
    for (i.run in 0:(runs-1)) {
      contents <- c(contents, paste(dps.path,
                                    args.string,
                                    sprintf("--pconj %.4e %s",
                                            pconj.values[i.pconj+1],
                                            sprintf("/home/kkentzo/%s/results.%d.%d.h5 &\n",
                                                    label, i.pconj, i.run)),
                                    sep=" "))

      job.id <- job.id + 1

      ## insert a wait??
      if (job.id %% cores == 0) {
        contents <- c(contents, "wait\n")
        wait.inserted <- T
      } else {
        wait.inserted <- F
      }
    }

  ## output number of jobs
  print(sprintf("Created %d jobs", job.id))

  ## write footer
  if (wait.inserted)
    contents <- c(contents, "\ndate\n")
  else
    contents <- c(contents, "wait\n\ndate\n")

  ## write contents to file
  if (is.na(fname))
    fname <- paste(label, "sh", sep=".")

  write(contents, fname, sep="\n")
  
}






## create the ht.global experiment
dps.create.experiment.script <- function(fname="htg.sh") {

  dps.create.cluster.script(label="htg", fname=fname, cores=200, runs=100,
                            ## == use these values for the extended simulations ==
                            ## pconj.values=c(10 ^ seq(-4, -1, length.out=30),
                            ##   seq(0.2, 0.5, 0.1)),
                            pconj.values=10 ^ seq(-4, -1, length.out=30),
                            beta=0.4, kappa=1, alpha=1, mu=5e-3,
                            mutate="bka", steps=250000, fparts=20)
  
}




## create a script pp.sh in the CWD that will post-process cluster results
## place it to the cluster directory that contains the results and
## pass its name to the sbatch command
dps.create.cluster.pp.script <- function( cores ) {

  contents <- paste("#!/bin/sh\n",
                    "#SBATCH --job-name=PP",
                    "#SBATCH --output=pp.out",
                    "#SBATCH --error=pp.err",
                    sprintf("#SBATCH --ntasks=%d\n", cores),
                    "date\n",
                    paste("LD_LIBRARY_PATH=/home/kkentzo/local/lib Rscript",
                          sprintf("/home/kkentzo/Projects/programs/dps.r . %d", cores)),
                    "\nwait\n\ndate\n",
                    sep="\n")


  write(contents, "pp.sh")
  
}






dps.create.cluster.plot.script <- function( n.pconj, n.runs ) {

  cores <- n.pconj * n.runs
  cmd <- paste("LD_LIBRARY_PATH=/home/kkentzo/local/lib",
               "Rscript",
               "/home/kkentzo/Projects/programs/dps.r plot")

  contents <- paste("#!/bin/sh\n",
                    "#SBATCH --job-name=BK.PLOT",
                    "#SBATCH --output=plot.out",
                    "#SBATCH --error=plot.err",
                    sprintf("#SBATCH --ntasks=%d\n", cores),
                    "date\n",
                    sep="\n")
  
  for (i.pconj in seq(0,n.pconj-1))
    for (i.run in seq(0,n.runs-1))
      contents <- paste(contents, sprintf("%s results.%d.%d.h5 figures 1 &\n",
                                          cmd, i.pconj, i.run),
                        sep="\n")

  contents <- paste(contents, "\nwait\n\ndate\n", sep="\n")

  ##print(contents)

  write(contents, "plot.new.sh")

  
}



## ==================================================================================
## ==================================================================================
##                              AD-HOC FUNCTIONS
## ==================================================================================
## ==================================================================================



ttest <- function() {

  alpha <- seq(-0.5,0.5,0.01)
  kappa <- seq(0,1,0.01)

  beta <- 0.4

  omega <- 1.5

  r <- matrix(NA, nrow=length(alpha), ncol=length(kappa))

  for (i in 1:length(alpha))
    for (j in 1:length(kappa))
      r[i,j] <- beta / (1. + kappa[j] * (3.5 + alpha[i]) / omega)

  contour(kappa, alpha, r, xlab=expression(kappa), ylab=expression(alpha))

  
}



tttest <- function() {

  n <- 1000

  kappa <- 0.3
  alpha <- 0.6
  
  delta.alpha <- rnorm(n, 0, 0.05)
  delta.kappa <- rnorm(n, 0, 0.05)

  omega <- 1.5

  f <- NULL

  for (i in 1:n)

    f <- c(f, 1. / (1. + ((kappa + delta.kappa[i]) * (7 * alpha + delta.alpha[i]) / omega)))


  plot(delta.alpha, f)

  print(cor(delta.alpha, f))
  
}




dps.pp.compare <- function(results.b, results.bka) {

  layout(matrix(c(1:4, rep(5,2), 6:7), byrow=T, ncol=2))

  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$div.inf - results.b$M$custom$death,
              results.bka$M$custom$div.inf - results.bka$M$custom$death),
        main="Host Growth", xlab=expression(p[c]), ylab="",
        log.take="x", col=c("blue", "red"))
  legend("left", c("NO-CNC", "CNC"),
         lwd=1, col=c("blue", "red"))

  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$death,
              results.bka$M$custom$death),
        main="Host Death", xlab=expression(p[c]), ylab="",
        log.take="x", col=c("blue", "red"))

  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$cn, results.bka$M$custom$cn),
        main="Copy Number", xlab=expression(p[c]), ylab="",
        log.take="x", col=c("blue", "red"))

  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$inf, results.bka$M$custom$inf),
        main="Infection", xlab=expression(p[c]), ylab="",
        log.take="x", col=c("blue", "red"))

  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$loss.cn, results.bka$M$custom$loss.cn),
        main="Seg Loss (per plasmid)", xlab=expression(p[c]), ylab="",
        log.take="xy", col=c("blue", "red"))

  mplot(cbind(results.b$pconj, results.bka$pconj),
        cbind(results.b$M$custom$rep.cn, results.bka$M$custom$rep.cn),
        main="Replication Rate (per plasmid)", ylab="", xlab=expression(p[c]),
        log.take="x")
  
  mplot(cbind(results.b$pconj, results.bka$pconj),
        ##cbind(results.b$M$global$M$ht, results.bka$M$global$M$ht),
        cbind(results.b$M$custom$ht.cn, results.bka$M$custom$ht.cn),
        main="Conjugation Rate (per plasmid)", ylab="", xlab=expression(p[c]),
        log.take="x")

  
  layout(matrix(1))
  
}




## compares OO relatedness between NO-CNC and CNC simulations
dps.pp.compare.relatedness <- function(results.b, results.bka) {

  var.names <- c("beta", "kappa", "alpha")

  ## form data frames
  r.nocnc.mean <- data.frame(Reduce(cbind, results.b$M$relatedness$oo, c()))
  r.cnc.mean <- data.frame(Reduce(cbind, results.bka$M$relatedness$oo, c()))

  ## rename columns
  names(r.nocnc.mean) <- var.names
  names(r.cnc.mean) <- var.names

  layout(matrix(1:3, nrow=1))

  for (var.name in var.names) {

    mplot(cbind(results.b$pconj, results.bka$pconj),
          cbind(r.nocnc.mean[[var.name]], r.cnc.mean[[var.name]]),
          xlab=expression(p[c]), ylab="", type="l",
          log.take="x", main=sprintf("Relatedness (%s)", var.name),
          col=c("blue", "red"))

    if (var.name == "beta") {
      legend("bottomleft", c("NO-CNC", "CNC"),
             lwd=1, col=c("blue", "red"))
    }
  }

  layout(matrix(1))
  
}



dps.pp.plot.price.beta <- function(results.b) {

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
           c("Total Selection (intra)",
             expression(paste("Selection due to ", n^R, " (intra)")),
             expression(paste("Selection due to ", n^H, " (intra)")),
             expression(paste("Selection due to ", n^D, " (intra)"))),
           lwd=1, lty=1:4, bty="n", cex=2, col="red")
    legend("bottomleft",
           c("Total Selection (inter)",
             expression(paste("Selection due to ", n^R, " (inter)")),
             expression(paste("Selection due to ", n^H, " (inter)")),
             expression(paste("Selection due to ", n^D, " (inter)"))),
           lwd=1, lty=1:4, bty="n", cex=2, col="blue")
    legend("topright", "Transmission Bias", ##inset=c(0, 0.35),
           lwd=1, bty="n", cex=2, col="green")
    
  
  
}



process.relatedness.adhoc <- function(path="relatedness") {

  for (fname in c("0.12", "1.13", "2.8")) {

    full.fname <- file.path(path, sprintf("results.%s.h5", fname))

    r <- dps.load(full.fname)

    dz <- dps.calc.price(r, 5000)

    save(dz, file=file.path(path, sprintf("dz.%s.xdr", fname)))
    
  }
  
}



## ============== MAIN SCRIPT ====================


## OPTIONS
## use:
##        Rscript dps.r PATH CORES
## to post-process the results of an experiment in PATH and produce results.xdr
##        Rscript dps.r plot RESULTS_PATH PNG_PATH CORES
## runs dps.analyze for all h5 files in RESULTS_PATH and saves the generated
## png files in PNG_PATH


## grab options
options <- commandArgs(trailingOnly=T)

if (length(options) == 2) {

  ## post-process the results
  dps.pp.experiment.parallel(options[1], as.numeric(options[2]))

} else if (length(options) == 4 && options[1] == "plot") {

  ## option[2] is the path to the h5 files [or FNAME for the serial version of plot.all()]
  ## option[3] is the path to the generated images
  ## option[4] is the number of cores to use
  dps.plot.all(options[2], options[3], as.numeric(options[4]))
  
}


