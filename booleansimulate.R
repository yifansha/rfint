library(iRF)

runSim <- function(n, p, x.params, y.params, dir, seed, n.cores=1) {
  set.seed(seed)
  x <- genX(n, x.params)
  y <- genY(x, y.params)
  ns <- ifelse(is.factor(y), 1, 10)
  if (!is.null(x.params$cov)) x.params$cov <- NULL
  
  int.return <- 1:5
  fit <- iRF(y=y, x=x, 
             interactions.return=int.return,
             n.iter=max(int.return),
             nodesize=ns,
             n.bootstrap=20,
             ntree=500,
             n.core=n.cores)
  
  gini.x <- decreaseGiniWrapper(x, y)
  file <- paste0(dir, 'iter', seed, '.Rdata')
  save(file=file, fit, gini.x, x.params, y.params)
}

genX <- function(n, x.params) {
  p <- x.params$p
  x <- matrix(rcauchy(n*p), nrow=n)
  return(x)
}

genOR <- function(x, tt, active.set=1:4) {
  stopifnot(ncol(x) >= max(active.set))
  y <- apply(x[,active.set] > tt, MAR=1, any)
  return(y)
}

genY <- function(x, y.params) {
  n <- nrow(x)
  p <- ncol(x)
  
  y <- genOR(x, y.params$thresh, y.params$active)
  
  if (y.params$noise == 0) {
    y <- as.factor(as.numeric(y))
  } else if (y.params$classification) {
    n.class <- trunc(table(y) * y.params$noise)
    idcs.swap <- mapply(function(cc, nn)
      sample(which(y == cc), nn, replace=TRUE),
      as.factor(names(n.class)), n.class)
    idcs.swap <- unlist(idcs.swap)
    y[idcs.swap] <- ifelse(y[idcs.swap], FALSE, TRUE)
    y <- as.factor(as.numeric(y))
  } else {
    y <- as.numeric(y) + rnorm(n, sd=y.params$noise)
  }
  return(y)
}


gini <- function(y) {
  if (is.factor(y)) yy <- as.numeric(y) - 1
  else yy <- y
  return(1 - mean(yy) ^ 2 - mean(1 - yy) ^ 2)
}

decreaseGini <- function(x, y) {
  n <- length(y)
  gini.y <- gini(y)
  x.order <- order(x)
  y.sort <- y[x.order]
  split.gini <- sapply(1:(n-1), function(i) 
    c(gini(y.sort[1:i]), gini(y.sort[-(1:i)])))
  weights <- rbind(1:(n-1), (n-1):1)
  child.gini <- colSums(weights * split.gini) / n
  dec.gini <- gini.y - min(child.gini)
  
  return(c(y.gini=gini.y, child.gini=min(child.gini), dec.gini=dec.gini))
}

decreaseGiniWrapper <- function(x, y) {
  gg <- apply(x, MAR=2, decreaseGini, y)
  return(gg)
}


n.cores <- 24

x.params <- list()
x.params$dist <- 'cauchy'
y.params <- list()
y.params$classification <- TRUE

thresh <- 3.2
noises <- c(0.1, 0.15, 0.2)
p <- 50
n <- 100
active <- 1:5


for (noise in noises) {
  print(paste('P:', p, 'N:', n))
  x.params$p <- p
  y.params$noise <- noise
  y.params$thresh <- thresh
  
  y.params$thresh <- thresh
  rr <- paste0('or', '_', y.params$thresh, '_Identity')
  y.params$active.set <- active
  y.params$k <- length(active)

  dir <- paste0('./booleanSims/noise_', noise, '/', rr, '/', 'nObs_', n, '/')
  dir.create(dir, recursive=TRUE)
  iters <- 1:2
  lapply(iters, function(ii) runSim(n, p, x.params, y.params, dir, ii, n.cores))
}






