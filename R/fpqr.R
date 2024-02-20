fpqr <- function(y, x, h, tau, nby, nbx, gpy, gpx, qc.type = c("dodge","choi","li")){


  qc.type <- match.arg(qc.type)

  BS.sol.x <- getAmat(data = x, nbf = nbx, gp = gpx)
  x1 <- BS.sol.x$Amat
  BS.sol.y <- getAmat(data = y, nbf = nby, gp = gpy)
  y1 <- BS.sol.y$Amat

  m.pqr <- pqr(y = y1, x = x1, tau = tau, h = h, qc.type = qc.type)
  fits <- (cbind(1,x1) %*% m.pqr$d.coef) %*% solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)

  b.hat <- BS.sol.x$evalbase %*% solve(BS.sol.x$sinp_mat) %*% m.pqr$d.coef[-1,] %*%
    solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)
  b0.hat <- t(solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)) %*% as.matrix(m.pqr$d.coef[1,])

  model.details <- list()
  model.details$nbx <- nbx
  model.details$gpx <- gpx
  model.details$m.pqr <- m.pqr
  model.details$BS.sol.x <- BS.sol.x
  model.details$BS.sol.y <- BS.sol.y

  return(list(fitted.values = fits, b0.hat = b0.hat, b.hat = b.hat, mdts = model.details))

}
