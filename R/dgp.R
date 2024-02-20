dgp <- function(n, nphi = 10, gpy, gpx, sd.error){
  
  X.phi1 <- function(nphi, gpx) {
    phi <- matrix(0, nphi, length(gpx))
    for (j in 1 : nphi)
      phi[j, ] <- (j^{-2}) * sqrt(2) * sin(j * pi * gpx)
    return(phi)
  }
  
  X.phi2 <- function(nphi, gpx) {
    phi <- matrix(0, nphi, length(gpx))
    for (j in 1 : nphi)
      phi[j, ] <- (j^{-2}) * sqrt(2) * cos(j * pi * gpx)
    return(phi)
  }
  
  rX.s <- function(nphi, gpx){
    xsi <- rnorm(2 * nphi, 0, 1)
    X <- xsi * rbind(X.phi1(nphi, gpx), X.phi2(nphi, gpx)) 
    Xs <- colSums(X)
    return(Xs)
  }
  
  alpha <- function(t) 2 * exp(-(t - 1)^2)
  beta <- function(t, s) {
    4 * cos(2 * pi * t) * sin(pi * s)
  }
  
  data <- list()
  coeffun <- list()
  s <- gpx
  t <- gpy
  ngpx <- length(gpx)
  ngpy <- length(gpy)
  alpha.t <- matrix(alpha(t), nrow = n, ncol = ngpy, byrow = TRUE)
  beta.ts <- outer(t, s, beta)
  data$X <- I(t(replicate(n, rX.s(nphi, s))))
  Xbeta.t <- (data$X) %*% t(beta.ts) * (s[2] - s[1])
  data$Xbeta <- I(Xbeta.t)
  eps <- matrix(rnorm(n * ngpy, 0, sd.error), n, ngpy)
  data$Ytrue <- I(alpha.t + Xbeta.t)
  data$Y <- I(alpha.t + Xbeta.t + eps)
  
  x <- matrix(data$X, ncol = ncol(data$X))
  y <- matrix(data$Y, ncol = ncol(data$Y))
  y.true <- matrix(data$Ytrue, ncol = ncol(data$Ytrue))
  
  coef.a <- alpha(gpy)
  coef.b <- outer(gpy, gpx, beta)
  return(list(x=x, y=y, y.true = y.true, coef.a = coef.a, coef.b = coef.b))
}
