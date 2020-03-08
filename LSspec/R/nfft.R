# nfft.R
# R code wrapper for nfft3 for fast computation of the Lomb-Scargle periodogram
# Aaron Springford
#####################

nfft <- function( t, x, N = 2*length(t), os = 1, setN = FALSE, tRange = NULL ){
  # t is a vector of times
  # x is a vector of values
  # N is the number of fourier coefs
  # os is the oversampling ratio (analogous to zero-pad ratio)
  # setN is whether to set N based on os value
  # tRange is the range of possible times. If NULL, the times are left as-is.
  #     If not NULL, an adjustment is made when scaling the times to [-0.5,0.5)
  
  if( length(t) != length(x) ) stop( "Arguments t and x are not the same length. \n" )
  m <- length(t)
  
  if( setN ) N <- 2*m*os
  
  if( N%%2 != 0 ) stop( "The number of fourier coefficients N must be even. \n")
  
  # time values need to be in [-0.5,0.5)
  nR <- c(-0.5 + 1e-5, - 0.5 + (1 - 1e-5)/os )
  # nR <- c(-0.5, -0.5 + (1 - 1/m)/os )
  if( is.null(tRange) ){
    tScaled <- scaleTimes(t, newRange = nR)
  }else{
    tScaled <- scaleTimes(t, newRange = nR, tRange)
  }
  
  # the nfft routine takes arguments
  #  t are the times
  #  y are the time-domain values
  #  N is the number of Fourier coefs
  #  M is the number of non-equispaced nodes (times)
  #  d are the frequency-domain values
  
  fftout <- .C( "nfft", as.double(tScaled), as.double(x), as.integer(N), as.integer(m), complex(N) )[[5]]
  fftout <- Conj(fftout[ -(1:(length(fftout)/2)) ])
  odds <- (1:length(fftout))%%2 == 0  # Because in R we are indexing from 1, not 0 
  fftout[ odds ] <- -fftout[ odds ]

fftout
}

scaleTimes <- function( x, newRange = c(-0.5, 0.5 - 1/(length(x))), r.x = range(x) ){
  sp.x <- r.x[2] - r.x[1]
(x-r.x[1])/sp.x*(newRange[2]-newRange[1]) + newRange[1]  
}


nfft2 <- function( t, x, N = 2*length(t) ){
  # t is a vector of times
  # x is a vector of values
  # N is the number of fourier coefs
  
  if( length(t) != length(x) ) stop( "Arguments t and x are not the same length. \n" )
  m <- length(t)
  
  if( N%%2 != 0 ) stop( "The number of fourier coefficients N must be even. \n")
  
  # time values need to be in [-0.5,0.5)
  tR = diff( range( t ) )
  a <- 0.5 - 1e-5
  tScaled <- 2*a*(t - min(t)) / tR - a
  
  # the nfft routine takes arguments
  #  t are the times
  #  y are the time-domain values
  #  N is the number of Fourier coefs
  #  M is the number of non-equispaced nodes (times)
  #  d are the frequency-domain values
  
  fftout <- .C( "nfft", as.double(tScaled), as.double(x), as.integer(N), as.integer(m), complex(N) )[[5]]
  #c( fftout[ -(1:(N/2)) ], Conj( fftout[1] ) )
  browser()
  fftout
}
