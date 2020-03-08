# ls.R
# computing the lomb-scargle periodogram
# using the nfft

PspecMT <- function( t, x, os = 1, tRange = range(t), w, k, tpI = NULL, subtract.mean = TRUE ){
  # Multitapering of the Schuster periodogram
  if( is.null(tpI) ) tpI <- dpssApproxFun( t, w, k, tRange = tRange )
  if( subtract.mean ) x <- x - mean(x)
  T <- diff( tRange ); N <- length(t)
  
  # The interpolated tapers
  tapers <- lapply( tpI, do.call, list(t) )
  # The tapers * data
  pf <- function( taper, t, x, T, dt ){
    taper/sum(taper^2) * x * sqrt( T/dt ) 
  }
  tpx <- lapply( tapers, pf, t = t, x = x, T = T, dt = T/N )
  
  # Compute nfft
  nfft.2 <- function( x, t, os, setN, tRange ) nfft( t, x, os = os, setN = TRUE, tRange = tRange )
  spc <- lapply( tpx, nfft.2, t = t, os = os, tRange = tRange )
  freq <- seq(from = 0, by = 1/os/diff(tRange), along.with = spc[[1]] )
  Pall <- sapply( spc, function(x) abs(x)^2 )
  
  # Averaging to get multitaper statistic
  # Note that we multiply by deltat * oversample
  out <- data.frame( freq = freq, P = rowMeans(Pall)* 0.5 * T / N / N )
  attr( out, "w" ) <- w
  attr( out, "k" ) <- k
  attr( out, "n" ) <- N
  attr( out, "tRange" ) <- tRange
  attr( out, "os" ) <- os
  attr( out, "subtract.mean" ) <- subtract.mean
  out
}



LSspec <- function( x, t, os = 1, tRange = range(t) ){

  # compute nfft of irregular time series
  sp <- nfft( t, x, os = os, setN = TRUE, tRange = tRange )
  # compute nfft of 1's at frequencies 2w
  spf <- nfft( t, rep(1, length(t)), N = 4*length(t)*os, os = os, tRange = tRange )
  spf <- spf[ rep_len(c(FALSE,TRUE),length(spf)) ]
  
  # compute components of LS periodogram from Press & Rybicki
  Sh <- -Im(sp)    # \sum_i x_i sin(2\pi f t_i)
  Ch <- Re(sp)     # \sum_i x_i cos(2\pi f t_i)
  S2 <- -Im(spf)   # \sum_i 1 sin(2\pi 2f)
  C2 <- Re(spf)    # \sum_i 1 cos(2\pi 2f)
  
  cos2wt <- C2/Mod(spf)
  sin2wt <- S2/Mod(spf)
  coswt <- sign(sin2wt)*sqrt(0.5*(1+cos2wt))
  sinwt <- sqrt(0.5*(1-cos2wt))
  
  
  A1 <- Ch*coswt + Sh*sinwt
  B1 <- Sh*coswt - Ch*sinwt
  A2 <- length(t)/2 + 0.5*C2*cos2wt + 0.5*S2*sin2wt
  B2 <- length(t) - A2
  
  P <- 0.5 * ( A1*A1/A2 + B1*B1/B2 )
  P <- P * 0.5 * diff(tRange) / length(t)
  freq <- seq(from = 0, by = 1/os/diff(tRange), along.with = P )
#browser()
# list( P, Sh, Ch, S2, C2, cos2wt, sin2wt, coswt, sinwt, A1, B1, A2, B2 )  
#data.frame( freq = seq(from = 0, by = 1/os/length(t), along.with = P ), P = P )
data.frame( freq = freq, P = P )
}

dpssApproxFun <- function( t, w, k, tRange = range(t) ){
  # t are the times
  # w is the bandwith
  # k is the number of tapers to return
  # tRange is the total range of times (possibly beyond range(t))
  require( multitaper )
  
  n <- length(t)
  T <- diff( tRange )
  # sample the dpss, then interpolate to needed values
  # The dpss is sampled regularly for accuracy (tridiagonal and whatnot)
  tp <- dpss( n = n, k = k, nw = n*w, returnEigenvalues = FALSE )$v
  tp <- as.data.frame(tp)
  # Now a cubic spline is fit to the regularly spaced dpss
  # splF <- function( taper, t ) smooth.spline( x = t, y = taper, all.knots = TRUE, keep.data = FALSE )
  splF <- function( taper, t ) splinefun( t, taper )

  out <- lapply( tp, splF, t = seq(from=tRange[1], to=tRange[2], length.out = n) )
  attr( out, "w" ) <- w
  attr( out, "k" ) <- k
  attr( out, "n" ) <- n
  attr( out, "tRange" ) <- tRange
out
}



LSspecMT <- function( t, x, os = 1, tRange = range(t), w, k, tpI = NULL, subtract.mean = TRUE ){
  # Multitapering of the Lomb-Scargle periodogram
  if( is.null(tpI) ) tpI <- dpssApproxFun( t, w, k, tRange = tRange )
  if( subtract.mean ) x <- x - mean(x)
  T <- diff( tRange ); N <- length(t)
  
  # The interpolated tapers * data
  # pf <- function( obj, t, x, T, dt ){
  #   taper <- predict( obj, t )$y
  #   # renormalize
  #   taper/sum(taper^2) * (x - mean(x)) * sqrt( T/dt ) 
  # }
  # tpx <- lapply( tpI, pf, t = t, x = x, T = T, dt = T/N )

  # The interpolated tapers
  tapers <- lapply( tpI, do.call, list(t) )
  # The tapers * data
  pf <- function( taper, t, x, T, dt ){
    taper/sum(taper^2) * x * sqrt( T/dt ) 
  }
  tpx <- lapply( tapers, pf, t = t, x = x, T = T, dt = T/N )
  
  # Compute spectra
  spc <- lapply( tpx, LSspec, t = t, os = os, tRange = tRange )
  freq <- spc[[1]]$freq
  Pall <- sapply( spc, getElement, name = "P" )
  
  # Averaging to get multitaper statistic
  # Note that we multiply by deltat * oversample
  out <- data.frame( freq = freq, P = rowMeans(Pall) )
  attr( out, "w" ) <- w
  attr( out, "k" ) <- k
  attr( out, "n" ) <- N
  attr( out, "tRange" ) <- tRange
  attr( out, "os" ) <- os
  attr( out, "subtract.mean" ) <- subtract.mean
out
}


