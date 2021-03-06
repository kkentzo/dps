suppressMessages(library(tseries))
suppressMessages(library(hdf5))
suppressMessages(library(gplots))
suppressMessages(library(multicore))

trySource <- function(path) {suppressWarnings(source(path)); cat("Sourcing", path, "\n");}

## insert contents of krutils.r
tryCatch(trySource('./krutils.r'),
         error=function(e) {
           tryCatch(trySource('R/krutils.r'),
                    error=function(e) cat("Unable to source krutils\n"))
         })

## *** See bottom of file for command line options ***




## =================================================================
##                 LOAD/SAVE/PROCESS HDF FILES
## =================================================================



## load the data from file FNAME
dps.load <- function( fname, verbose=F ) {

  if (verbose)
    print(sprintf("Loading file %s", fname))

  results <- hdf5load(fname, load=F, tidy=T)

  results
  

}




## =================================================================
## save the results to file FNAME
dps.save <- function( results, fname ) {

  dynamics <- results$dynamics
  settings <- results$settings
  population <- results$population
  histograms <- results$histograms

  hdf5save(fname, "dynamics", "settings", "population", "histograms")
  
}





## =================================================================
##                     PLOT RESULTS (DYNAMICS/HISTS)
## =================================================================



## =================================================================
## calculates and returns the association between the pair of
## variables specified by NAMES
## NORM.BY will divide the result by the variance of none, one 
## or both components (for "xy" the correlation coefficient is computed)
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


## =================================================================
## loads and returns the saved dz object calculated by dps.calc.price()
dps.load.price <- function(fname) {

  load(fname)
  dz
  
}


## saves the DZ object calculated by dps.calc.price()
dps.save.price <- function(dz, fname) {

  save(dz, file=fname, compress=T)
  
}




## =================================================================
## calculates and returns the 3 components of the price equation
## for beta, kappa and alpha (uses 3 cores to do that in parallel)
## STEPS.RANGE is a 2-tuple
dps.calc.price <- function(results, window.size=500, steps.range=NA, 
                           decompose=T, relatedness=T,
                           variance=T, correlations=F)
{

  if (length(steps.range) != 2) 
    take.steps <- 1:length(results$dynamics$global$M$fitness)
  else
    take.steps <- steps.range[1]:steps.range[2]
  
  dz <- list(fitness=ma(results$dynamics$global$M$fitness[take.steps],
                 window.size=window.size),
             host.growth=ma(results$dynamics$counters$div.all /
                 results$dynamics$counters$n, window.size),
             window.size=window.size)

  ## returns a list with members: z.name, and if decompose: "nr", "ht", "death"
  calc.components.for.level <- function(z.name, level) {
    l <- list()
    l[[level]] <- ma(dps.calc.association(results$dynamics[[level]],
                                          paste(z.name, "fitness", sep="."),
                                          norm.by="",
                                          steps.range=steps.range),
                     window.size=window.size) / dz$fitness
    if (decompose) {
      for (name in c("nr", "ht", "death")) {
        l[[paste(level, name, sep=".")]] <-
          ma(dps.calc.association(results$dynamics[[level]],
                                  paste(z.name, name, sep="."),
                                  norm.by="",
                                  steps.range=steps.range),
             window.size=window.size) / dz$fitness
      }
    }
    l
  }

  ## returns a list with all the price components calculated
  ## for variable Z.NAME
  calc.price.for.var <- function(z.name) {
    ## calculate and store transmission bias
    r <- data.frame(tbias=ma(
                        results$dynamics$global$M[[sprintf("t%s", z.name)]][take.steps],
                        window.size=window.size) / dz$fitness)

    ## calculate covariances
    level.names <- c("global", "inter", "intra")
    vars <- mclapply(level.names,
                     function(level) calc.components.for.level(z.name, level),
                     mc.preschedule=F, mc.cores=length(level.names))
    ## attach all the results to data.frame r
    for (v in vars)
      for (v.name in names(v))
        r[[v.name]] <- v[[v.name]]
    ## Calculate total selection
    r$total <- r$global + r$tbias
    ## Calculate relatedness (separate components for cov and var)
    if (relatedness) {
      r$r <- data.frame(
          cov=ma(results$dynamics$relatedness$oo$cov[[z.name]][take.steps],
              window.size=window.size),
          var=ma(results$dynamics$relatedness$oo$var[[z.name]][take.steps],
              window.size=window.size))
      r$r.wg <- data.frame(
          cov=ma(results$dynamics$relatedness$wg$cov[[z.name]][take.steps],
              window.size=window.size),
          var=ma(results$dynamics$relatedness$wg$var[[z.name]][take.steps],
              window.size=window.size))
    }

    ## return r
    r
    
  }

  ## calculate price equation components
  z.names <- c("beta", "kappa", "alpha")
  vars <- mclapply(z.names, calc.price.for.var,
                   mc.preschedule=F, mc.cores=length(z.names))
  names(vars) <- z.names
  dz <- c(dz, vars)

  ## calculate correlations?
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

  ## calculate moving averages for variances?
  if (variance)
    dz$variance <- dps.calc.variance(results, window.size, steps.range)
  
  dz
  
}



## =========================================================================
## calculate the moving averages of the plasmid parameters variances
## returns a list with the data.frames global, inter, intra
dps.calc.variance <- function(results, window.size=500, steps.range=NA) {

  if (length(steps.range) != 2) 
    take.steps <- 1:length(results$dynamics$global$M$fitness)
  else
    take.steps <- steps.range[1]:steps.range[2]

  level.names <- c("global", "inter", "intra")

  l <-
    mclapply(level.names,
             function(level.name) {
               z.names <- c("beta", "kappa", "alpha")
               ll <- mclapply(z.names,
                              function(z.name) {
                                ma(results$dynamics[[level.name]]$V[[z.name]][take.steps],
                                   window.size=window.size)
                              }, mc.cores=length(z.names))
               names(ll) <- z.names
               as.data.frame(ll)
             }, mc.cores=length(level.names))

  names(l) <- level.names

  l
  
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
## plot the dynamics of price equation events components
dps.plot.price.events <- function(dz, events=c("nr", "death"),
                                  z.names=c("beta", "kappa", "alpha"),
                                  steps.range=NA) {

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dz$beta)
  else
    steps.range <- steps.range[1]:steps.range[2]

  level.names <- c("global", "inter", "intra")
  layout(matrix(1:(3*length(z.names)), ncol=3, byrow=T))

  ## plot components of selection due to the specified events
  for (z.name in z.names) {
    for (level in level.names) {
      ## plot time series
      nr <- dz[[z.name]][[paste(level,"nr",sep=".")]][steps.range]
      nd <- - dz[[z.name]][[paste(level,"death",sep=".")]][steps.range]
      total <- nr + nd
      mplot(steps.range, cbind(total, nr, nd), xlab="Time",
            main=paste(z.name, level, sep=" "),
            col=c("black", "blue", "red"))
      abline(h=0)
      legend("bottomleft", c("total", "nr", "death"),
             lwd=1, col=c("black", "blue", "red"))
    }
  }

  layout(matrix(1))
  
}




## =========================================================================
## plot the dynamics of price equation events components
dps.plot.price.events2 <- function(dz, events=c("nr", "death"),
                                   z.names=c("beta", "kappa", "alpha"),
                                   steps.range=NA) {

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dz$beta)
  else
    steps.range <- steps.range[1]:steps.range[2]

  layout(matrix(1:3, ncol=3, byrow=F))

  ## plot components of selection due to the specified events
  for (z.name in z.names) {
    ## plot dynamics
    ##mplot(steps.range, ma(results$dynamics$counters$rep[steps.range], 500))
    ## plot BHS
    nr <- dz[[z.name]]$inter.nr[steps.range]
    nd <- - dz[[z.name]]$inter.death[steps.range]
    total <- nr + nd
    mplot(steps.range, cbind(nr, nd), xlab="Time",
          main=sprintf("BHS on %s", z.name),
          col=c("black", "blue", "red"))
    abline(h=0)
    legend("bottomleft", c("total", "nr", "death"),
           lwd=1, col=c("black", "blue", "red"))
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
dps.plot.relatedness <- function( dz ) {

  layout(matrix(1:3, nrow=1))

  ## plot beta
  mplot(cbind(dz$beta$r$cov, dz$beta$r$var),
        xlab="Time", ylab="", main="beta", log.take="y")
  ## plot kappa
  mplot(cbind(dz$kappa$r$cov, dz$kappa$r$var),
        xlab="Time", ylab="", main="kappa", log.take="y")
  ## plot alpha
  mplot(cbind(dz$alpha$r$cov, dz$alpha$r$var),
        xlab="Time", ylab="", main="alpha", log.take="y")

  layout(matrix(1))

}


## =========================================================================
## plot the relatedness components (SORTED BY RELATEDNESS COMPONENTS)
dps.plot.relatedness2 <- function( dz ) {

  layout(matrix(1:3, nrow=1))

  ## plot relatedness
  mplot(cbind(dz$beta$r$cov / dz$beta$r$var,
              dz$kappa$r$cov / dz$kappa$r$var,
              dz$alpha$r$cov / dz$alpha$r$var),
        xlab="Evolutionary Time", ylab="", main="Relatedness")

  ## plot COV component
  mplot(cbind(dz$beta$r$cov,
              dz$kappa$r$cov,
              dz$alpha$r$cov),
        xlab="Evolutionary Time", ylab="", main="Covariance", log.take="y")

  ## plot VAR component
  mplot(cbind(dz$beta$r$var,
              dz$kappa$r$var,
              dz$alpha$r$var),
        xlab="Evolutionary Time", ylab="", main="Variance", log.take="y")

  layout(matrix(1))

}



## =================================================================
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



## ==================================================================
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








dps.plot.histogram <- function( results, name="age", fpart=NULL ) {

  x <- results$histograms[[name]]$x

  if (is.null(fpart)) 
    y <- colSums(results$histograms[[name]]$y)
  else
    y <- results$histograms[[name]]$y[fpart,]

  x <- x[y!=0]
  y <- y[y!=0]

  plot(x, log10(y), t="b", xlab=name, ylab="Frequency", yaxt="n")
  decorate.log("y")
  
  
}



dps.plot.cn.hists <- function(results) {

  x <- 0.5 + results$histograms$cn$x
  y <- t(results$histograms$cn$y)

  mplot(log10(x[2:length(x)]), y[2:length(x),], col="black", xaxt='n')
  decorate.log("x")
  abline(v=log10(7))
  
}






## ==================================================================
## CMD LINE FUNCTION : calls dps.analyze() and dps.plot.price() for
## all files in H5.PATH and saves the generated png files in PNG.PATH
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

      ## plot dynamics
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

  ## ================================================================
  ## run the plotting in parallel
  pvec(fnames, plot.results, mc.cores=cores)

}







## ==================================================================
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








## ==================================================================
## CMD LINE FUNCTION : called from the command line
## if CORES > 1, it calls dps.plot.all.parallel()
## otherwise, it calls dps.plot.all.single()
dps.plot.all <- function( h5.path, png.path, cores ) {

  if (cores > 1)
    dps.plot.all.parallel(h5.path, png.path, cores)
  else
    dps.plot.all.single(h5.path, png.path, 1)
  
}




## ==================================================================
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







### =================================================================
### =================================================================
###                      CLUSTER EXPERIMENTS
### =================================================================
### =================================================================





## load the XDR file with the post-processed results
dps.pp.load <- function( fname ) {

  load(fname)
  results
  
}








## ==================================================================
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



  ## ================================================================
  ## process the results for the specified index of pconj (for all runs)
  ## return a RESULTS list
  process.pconj <- function( i.pconj ) {

    pconj.values <- NA

    results <- list(counters=list(), ## automatic
                    relatedness=list(), ## automatic
                    drelatedness=list(), ## automatic
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
                        inf=array(NA, runs)),
                    histograms=list())

    
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
      for (r.type in names(r$dynamics$relatedness)) {
        for (name in names(r$dynamics$relatedness[[r.type]]$cov)) {
          ## populate relatedness coefficients in
          ## relatedness[[rtype]][[beta|kappa|alpha]]
          if (is.null(results$relatedness[[r.type]][[name]])) 
            results$relatedness[[r.type]][[name]] <- array(NA, runs)
          
          results$relatedness[[r.type]][[name]][i.run] <-
            mean(r$dynamics$relatedness[[r.type]]$cov[[name]][seq.steps] /
                 r$dynamics$relatedness[[r.type]]$var[[name]][seq.steps],
                 na.rm=T)
        }
      }

      ## process DRELATEDNESS (relatedness decomposed into cov and var)
      for (r.type in names(r$dynamics$relatedness)) { ## ==> oo, wg
        for (s.type in c("cov", "var")) { ## cov, var
          for (name in names(r$dynamics$relatedness[[r.type]]$cov)) { ## ==> b,k,a
            if (is.null(results$drelatedness[[r.type]][[s.type]][[name]])) 
              results$drelatedness[[r.type]][[s.type]][[name]] <- array(NA, runs)

            results$drelatedness[[r.type]][[s.type]][[name]][i.run] <-
              mean(r$dynamics$relatedness[[r.type]][[s.type]][[name]][seq.steps],
                   na.rm=T)
          }
        }
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

      ## process HISTOGRAMS
      for (name in names(r$histograms)) {
        if (is.null(results$histograms[[name]]))
          results$histograms[[name]] <- r$histograms[[name]]
        else if ("z" %in% names(results$histograms[[name]]))
          results$histograms[[name]]$z <-
            results$histograms[[name]]$z + r$histograms[[name]]$z
        else
          results$histograms[[name]]$y <-
            results$histograms[[name]]$y + r$histograms[[name]]$y
      }
      

      ## remove results from memory
      rm(r)
      
    }

    ## form aggregated results
    results.agg <- list(M=list(counters=list(),relatedness=list(),
                            drelatedness=list(), custom=list()),
                        S=list(counters=list(), relatedness=list(),
                            drelatedness=list(), custom=list()),
                        histograms=results$histograms,
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
        
      } else if (level == "drelatedness") {

        for (r.type in names(results[[level]])) { ## ==> oo, wg
          for (s.type in c("cov", "var")) { ## cov, var
            if (is.null(results.agg$M[[level]][[r.type]][[s.type]])) {
              results.agg$M[[level]][[r.type]][[s.type]] <- list()
              results.agg$S[[level]][[r.type]][[s.type]] <- list()
            }
            
            for (name in names(results[[level]][[r.type]][[s.type]])) { ## ==> b,k,a
              results.agg$M[[level]][[r.type]][[s.type]][[name]] <-
                mean(results[[level]][[r.type]][[s.type]][[name]], na.rm=T)
              results.agg$S[[level]][[r.type]][[s.type]][[name]] <-
                sd(results[[level]][[r.type]][[s.type]][[name]], na.rm=T)
            }
          }
        }
      } else {

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
  



  ## ================================================================
  ## AD-HOC TESTING (COMMENT OUT WHEN DONE)
  ## len.pconj <- 1
  ## results <- process.pconj( 1 )
  ## return(results)
  ## ================================================================

  
  
  ## ================================================================
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

    histograms <- list(results[["1"]]$histograms)

    for (i.pconj in 2:len.pconj) {
      ## print value of i.pconj
      print(sprintf("i=%d", i.pconj))
      ## first, the PCONJ.VALUES
      results[["1"]]$pconj.values <- c(results[["1"]]$pconj.values,
                                       results[[as.character(i.pconj)]]$pconj.values)
      ## second, store the HISTOGRAMS
      histograms[[i.pconj]] <- results[[as.character(i.pconj)]]$histograms
      ## third, the LEVELS (counters, global, inter, intra, relatedness)
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
          
        } else if (level == "drelatedness") {

          for (r.type in names(results[["1"]]$M[[level]])) { ## ==> oo, wg
            for (s.type in c("cov", "var")) { ## cov, var
              for (name in names(results[["1"]]$M[[level]][[r.type]]$cov)) { ## ==> b,k,a
                results[["1"]]$M[[level]][[r.type]][[s.type]][[name]] <-
                  c(results[["1"]]$M[[level]][[r.type]][[s.type]][[name]],
                    results[[as.character(i.pconj)]]$M[[level]][[r.type]][[s.type]][[name]])
                results[["1"]]$S[[level]][[r.type]][[s.type]][[name]] <-
                  c(results[["1"]]$S[[level]][[r.type]][[s.type]][[name]],
                    results[[as.character(i.pconj)]]$S[[level]][[r.type]][[s.type]][[name]])
              }
            }
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
    ## replace the histograms in the first element
    ## with the aggregated list
    results[["1"]]$histograms <- histograms
  }
  
  ## form results
  results <- results[["1"]]

  ## store number of runs
  results$runs <- runs

  ## save as HDR 
  save(results, file=file.path(path, "results.xdr"), compress=T)

}









## ==================================================================
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




## ==================================================================
## ==================================================================
##                 POST-PROCESS COMP EXPERIMENTS
## ==================================================================
## ==================================================================

## post-process the results of the comp experiment in PATH
## and save them in $PATH/results.xdr
dps.pp.comp <- function(path, cores=1, stable="beta") {

  ## gather all results files
  fnames <- Sys.glob(file.path(path, "results.[0-9]*.h5"))
  ## read in value of pconj
  r <- dps.load(fnames[[1]])
  pconj <- r$settings$pconj
  rm(r)
  ## count the the number of kappa, alpha and runs values
  len.x <- 0
  len.y <- 0
  runs <- 0
  for (fname in fnames) {
    ## extract numerical tokens 
    tokens1 <- strsplit(fname, "results\\.")[[1]]
    tokens <- strsplit(tokens1[length(tokens1)], "\\.")[[1]]
    x <- as.numeric(tokens[1])
    y <- as.numeric(tokens[2])
    r <- as.numeric(tokens[3])
    if (y == 0 && r == 0)
      len.x <- len.x + 1
    if (x == 0 && r == 0)
      len.y <- len.y + 1
    if (x == 0 && y == 0)
      runs <- runs + 1
  }

  cat("X=", len.x, "\nY=", len.y,
      "\nRUNS=", runs, "\npconj=", pconj, "\n\n")

  ## process all runs for the specific combination
  ## of (i.x, i.y) <-- these are 1-based indices
  process.x.y <- function( i.x, i.y ) {
    result <- list(mut.wins=rep(0,runs),
                   both.exist=rep(0,runs),
                   steps=rep(0,runs))
    for (i.run in seq(runs)) {
      ## load file
      r <- dps.load(file.path(path, sprintf("results.%d.%d.%d.h5",
                                            i.x-1, i.y-1, i.run-1)))
      ## record values of kappa and alpha
      if (i.run == 1) {
        if (stable == "beta") {
          x <- r$competition$contenders$B$kappa
        } else {
          x <- r$competition$contenders$B$beta
        }
        y <- r$competition$contenders$B$alpha
      }
      ## figure out (and store) number of steps
      result$steps[i.run] <- nrow(r$competition$frequencies)
      ## does mutant win??
      if (r$competition$frequencies$A[result$steps[i.run]] == 0 &&
          r$competition$frequencies$B[result$steps[i.run]] > 0)
        result$mut.wins[i.run] <- 1
      ## do both exist?
      if (r$competition$frequencies$A[result$steps[i.run]] > 0 &&
          r$competition$frequencies$B[result$steps[i.run]] > 0)
        result$both.exist[i.run] <- 1
      rm(r)
    }
    ## calculate means and return results
    c(lapply(result, mean),
      i.x=i.x, i.y=i.y, x=x, y=y)
  }


  ## formulate the list of jobs to be calculated
  jobs <- list()
  ctr <- 1
  for (i in seq(len.x)) {
    for (j in seq(len.y)) {
      jobs[[ctr]] <- list(i.x=i, i.y=j)
      ctr <- ctr + 1
    }
  }
  ## schedule the jobs in parallel
  job.results <- mclapply(jobs,
                          function(job) process.x.y(job$i.x, job$i.y),
                          mc.cores=cores)

  ## aggregate results
  results <- list(pconj=pconj,
                  runs=runs,
                  x.values=rep(NA, len.x),
                  y.values=rep(NA, len.y),
                  mut.wins=matrix(NA, nrow=len.x, ncol=len.y),
                  both.exist=matrix(NA, nrow=len.x, ncol=len.y),
                  steps=matrix(NA, nrow=len.x, ncol=len.y))
  for (r in job.results) {
    results$x.values[r$i.x] <- r$x
    results$y.values[r$i.y] <- r$y
    results$mut.wins[r$i.x,r$i.y] <- r$mut.wins
    results$both.exist[r$i.x,r$i.y] <- r$both.exist
    results$steps[r$i.x,r$i.y] <- r$steps
  }
  ## attach experiment label to results
  if (stable == "beta")
    results$label <- "ka"
  else if (stable == "kappa")
    results$label <- "ba"

  save(results, file=file.path(path, "results.xdr"), compress=T)
}



## loads the post-processed comp results from FNAME
dps.pp.comp.load <- function(fname) {
  load(fname)
  results
}


## plots the results of the supplied comp experiment
## the experiment is derived from results$label
dps.pp.comp.plot <- function(results, metric="mut.wins") {

  if (results$label == "ba") {
    xlab <- expression(beta)
    xlim <- c(0.30,0.5)
    ylim <- c(0.80,1.0)
    annot <- c(0.4,0.9)
  } else if (results$label == "ka") {
    xlab <- expression(kappa)
    xlim <- c(0.80,1)
    ylim <- c(0.80,1)
    annot <- c(0.9,0.9)
  } else {
    stop(sprintf("Experiment %s not supported\n", results$label))
  }

  z.rng <- range(results[[metric]])
  if (z.rng[1] <= 1) {
    levels <- seq(0, 0.4, by=0.01)
    nlevels <- length(levels)
  } else {
    nlevels = 20
    levels = pretty(z.rng, nlevels)
  }
  
  filled.contour(results$x.values, results$y.values, results[[metric]],
                 main=results$label, xlab=xlab, ylab=expression(alpha),
                 levels=levels, nlevels=nlevels, xlim=xlim, ylim=ylim,
                 col=colorpanel(length(levels), "white", "red",),##"grey10"),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   points(annot[1],annot[2])
                   text(annot[1], annot[2], labels="WT", cex=1.5, pos=4)
                 } )
}



## ==================================================================
## ==================================================================
##                       AD-HOC FUNCTIONS
## ==================================================================
## ==================================================================



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






## compares the cov and var components OO relatedness
## between NO-CNC and CNC simulations
dps.pp.compare.relatedness.components <- function(results.b, results.bka,
                                                  only=c("beta")) {

  var.names <- c("beta", "kappa", "alpha")
  if (is.null(only))
    only <- var.names

  ## form data frames
  cov.nocnc.mean <- data.frame(Reduce(cbind, results.b$M$drelatedness$oo$cov, c()))
  var.nocnc.mean <- data.frame(Reduce(cbind, results.b$M$drelatedness$oo$var, c()))
  
  cov.cnc.mean <- data.frame(Reduce(cbind, results.bka$M$drelatedness$oo$cov, c()))
  var.cnc.mean <- data.frame(Reduce(cbind, results.bka$M$drelatedness$oo$var, c()))

  ## rename columns
  names(cov.nocnc.mean) <- var.names
  names(var.nocnc.mean) <- var.names
  names(cov.cnc.mean) <- var.names
  names(var.cnc.mean) <- var.names

  layout(matrix(1:length(only), nrow=1))

  for (var.name in var.names) {

    if (var.name %in% only) {

      mplot(results.b$pconj,
            cbind(cov.nocnc.mean[[var.name]], cov.cnc.mean[[var.name]],
                  var.nocnc.mean[[var.name]], var.cnc.mean[[var.name]]),
            xlab=expression(p[c]), ylab="", type="l", ltype=c(1,1,2,2),
            log.take="xy", main=sprintf("Relatedness (%s)", var.name),
            col=c("blue", "red", "blue", "red"))

      if (var.name == "beta") {
        legend("left", c("NO-CNC", "CNC"),
               lwd=1, col=c("blue", "red"))
        legend("right", c("cov", "var"),
               lwd=1, lty=c(1,2), col="black")
      }
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
## to post-process the results of an experiment in PATH
## and produce $PATH/results.xdr
##        Rscript dps.r comp PATH CORES STABLE
## to post-process the results of a COMP experiment in PATH
## and produce PATH/results.xdr
##        Rscript dps.r plot RESULTS_PATH PNG_PATH CORES
## runs dps.analyze for all h5 files in RESULTS_PATH and saves the generated
## png files in PNG_PATH


## grab all options
all.options <- commandArgs(trailingOnly=F)
## figure out whether any args have been passed
idx.args <- grep("--args", all.options) + 1

if (length(idx.args) > 0) {

  ## discover the directory of `this` script
  script.name <- sub("--file=", "",
                     all.options[grep("--file=", all.options)])
  this.path <- dirname(script.name)
  ## try to source krutils
  trySource(file.path(this.path, "krutils.r"))
  ## isolate the program options
  options <- all.options[idx.args:length(all.options)]
  ## parse program options and decide on course of action
  if (length(options) == 2) {

    ## post-process the results
    dps.pp.experiment.parallel(options[1], as.numeric(options[2]))

  } else if (length(options) == 4 && options[1] == "comp") {
    dps.pp.comp(options[2], as.numeric(options[3]), options[4])
  } else if (length(options) == 4 && options[1] == "plot") {

    ## option[2] is the path to the h5 files [or FNAME for the serial version of plot.all()]
    ## option[3] is the path to the generated images
    ## option[4] is the number of cores to use
    dps.plot.all(options[2], options[3], as.numeric(options[4]))
    
  }
}
