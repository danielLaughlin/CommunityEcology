minEntropy <- function (t2c = NULL, constraints = NULL, t2d,
                        obj = c("QH","Q", "H"), phi = 0.5, capd = FALSE, euclid = TRUE,
                        pars=c(rep(1/N,N)) ) 
{
  if (!is.null(t2c)) {
    if (!inherits(t2c, "matrix") | !is.numeric(t2c)) {
      stop("object \"t2c\" is not a numerical matrix")
    }
    obj <- match.arg(obj)
    if (ncol(t2c) >= nrow(t2c) - 1) {
      stop("There are more traits than species. Use fewer traits or more species.")
    }
    if (any(is.na(t2c))) {
      stop("Traits to constrain \"t2c\" contains NA")
    }
    if (is.null(colnames(t2c))) 
      warning("\"t2c\" colname vector (traits) is empty")
    if (is.null(rownames(t2c))) 
      warning("\"t2c\" rowname vector (species) is empty")
  }
  if (phi < 0 | phi > 1 | !is.numeric(phi)) {
    stop("phi must be between 0 and 1")
  }
  if (inherits(t2d, "matrix") & is.numeric(t2d)) {
    if (is.null(colnames(t2d))) 
      warning("\"t2d\" colname vector (traits) is empty")
    if (is.null(rownames(t2d))) 
      warning("\"t2d\" rowname vector (species) is empty")
    d <- stats::dist(t2d, upper = TRUE, diag = TRUE)
  }
  else if (inherits(t2d, "dist")) {
    d <- t2d
    if (is.null(names(t2d))) 
      warning("\"t2d\" name vector (species) is empty")
  }
  else {
    stop("object \"t2d\" must be a numerical trait matrix or a \"dist\" matrix")
  }
  if (any(is.na(d))) 
    stop("Distance matrix \"t2d\" contains NA")
  oldw <- getOption("warn")
  options(warn = -1)
  test <- !ade4::is.euclid(d)
  options(warn = oldw)
  if (euclid & test) 
    d <- as.matrix(ade4::quasieuclid(stats::as.dist(d)))
  d <- as.matrix(d)
  if (!is.null(t2c)) {
    if (nrow(d) != nrow(t2c)) {
      stop("objects \"t2d\" and \"t2c\" do not describe the same number of species")
    }
    if (!is.null(row.names(d)) & !is.null(row.names(t2c)) & 
        !identical(row.names(d), row.names(t2c))) {
      stop("\"t2d\" and \"t2c\" species names do not match")
    }
  }
  if (!is.null(constraints)) {
    if (is.null(t2c)) {
      stop("If \"constraints\" is specified, \"t2c\" must be specified")
    }
    if (length(constraints) != ncol(t2c)) {
      stop("Number of trait constraints number must be equal to the number of traits to constrain.")
    }
    if (!is.numeric(constraints) | any(is.na(constraints))) {
      stop("Trait constraint must all be numerical non-missing values.")
    }
    range_t2c <- apply(t2c, 2, range)
    if (any(range_t2c[1, ] > constraints | range_t2c[2, ] < 
            constraints)) {
      stop("Trait constraint must be within the range of available trait values in the species pool.")
    }
    if (!is.null(names(constraints)) & !is.null(colnames(t2c)) & 
        !identical(names(constraints), colnames(t2c))) {
      stop("\"constraints\" and \"t2c\" trait names don't match")
    }
    if (is.null(names(constraints))) 
      warning("\"constraints\" name vector is empty")
  }
  N <- nrow(d)
  mean.dist <- mean(d[upper.tri(d, diag = FALSE)])
  dcap.fun <- function(d) {
    d[d > mean.dist] <- mean.dist
    return(d)
  }
  if (capd) {
    d <- dcap.fun(d)
  }
  mean.t2d <- mean(t2d)
  tdist <- abs(t2d - mean.t2d)
  tdist[tdist == 0] <- 0.001
  qh <- function(p) {
    -(t(p) %*% (d/2) %*% p * phi + -(t(p) %*% log(p)) * (1 - 
                                                           phi))
  }
  q <- function(p) {
    -(t(p) %*% (d/2) %*% p)
  }
  h <- function(p) {
    (-(t(p) %*% log(p))) # minimize entropy!
  }
  if (!is.null(constraints)) {
    eqfun <- function(p) {
      z1 <- sum(p)
      cwms <- t(as.matrix(p)) %*% t2c
      return(c(z1, cwms))
    }
    all.constraints <- c(1, constraints)
  }
  else {
    eqfun <- function(p) {
      z1 <- sum(p)
      return(z1)
    }
    all.constraints <- 1
  }
  oldw <- getOption("warn")
  options(warn = -1)
  if (obj == "QH") {
    res <- Rsolnp::solnp(pars = pars, fun = qh, 
                         eqfun, eqB = all.constraints, LB = c(rep(0, N)), 
                         UB = c(rep(1, N)))
  }
  else if (obj == "Q") {
    res <- Rsolnp::solnp(pars = pars, fun = q, 
                         eqfun, eqB = all.constraints, LB = c(rep(0, N)), 
                         UB = c(rep(1, N)))
  }
  else if (obj == "H") {
    res <- Rsolnp::solnp(pars = pars, fun = h, 
                         eqfun, eqB = all.constraints, LB = c(rep(0, N)), 
                         UB = c(rep(1, N)))
  }
  result = list()
  result$prob <- as.matrix(res$pars)
  if (inherits(t2d, "matrix")) {
    rownames(result$prob) <- row.names(t2d)
    result$cwm_t2d <- as.matrix(res$pars %*% t2d)
  }
  else {
    rownames(result$prob) <- names(t2d)
  }
  if (!is.null(t2c)) {
    result$cwm_t2c <- as.matrix(res$pars %*% t2c)
  }
  result$H <- h(res$pars) * -1
  result$Q <- q(res$pars) * -1
  result$objval <- as.matrix(res$values)
  result$lagrange <- as.matrix(res$lagrange)
  result$hessian <- as.matrix(res$hessian)
  result$convergence <- res$convergence
  options(warn = oldw)
  return(result)
}
