# SMPM shared simulation engine | base R only
# All functional vectors are stored in weighted Euclidean coordinates.
smpm_config <- function(M = 201L, alpha = .4, beta = .3,
                        p = rep(1/3, 3), inverse_M = 5001L) {
  stopifnot(M >= 11, inverse_M >= 1001, length(p) == 3L,
            all(p > 0), abs(sum(p)-1) < 1e-10)
  grid <- seq(0, 1, length.out = M)
  w <- rep(1/(M-1), M); w[c(1,M)] <- w[c(1,M)]/2
  a <- cbind(c(1,1,0,0,0)/sqrt(2), c(0,0,1,-1,0)/sqrt(2))
  theta <- rbind(c(alpha,beta), c(-alpha,beta), c(0,-2*beta))
  H <- function(t) cbind(sqrt(2)/(2*pi)*sin(2*pi*t),
                          sqrt(2)/(4*pi)*sin(4*pi*t))
  U <- function(t) cbind(sqrt(2)*cos(2*pi*t), sqrt(2)*cos(4*pi*t))
  # Layout: component 1 grid, component 2 grid, etc.
  vectorize <- function(F) as.vector(t(F * rep(sqrt(w), each=nrow(F))))
  means <- lapply(1:3, function(g) {
    F <- matrix(rep(grid, each=5L), nrow=5L)
    for (r in 1:5) F[r,] <- grid + rowSums(sweep(H(grid), 2, theta[g,]*a[r,], '*'))
    vectorize(F)
  })
  means <- do.call(rbind, means)
  z <- colSums(means * p)
  V <- sweep(means, 2, z, '-')
  A <- sweep(V, 1, sqrt(p), '*')
  K <- tcrossprod(A)
  ee <- eigen(K, symmetric=TRUE)
  stopifnot(ee$values[2] > 0)
  lambdas <- ee$values[1:2]
  phi <- t(A) %*% ee$vectors[,1:2,drop=FALSE]
  phi <- sweep(phi, 2, sqrt(lambdas), '/')
  # Analytical directions apply only to equal group probabilities.
  if (max(abs(p-1/3)) < 1e-12) {
    phi_exact <- cbind(vectorize(outer(a[,1], sqrt(2)*sin(2*pi*grid))),
                       vectorize(outer(a[,2], sqrt(2)*sin(4*pi*grid))))
    phi_exact <- sweep(phi_exact, 2, sqrt(colSums(phi_exact^2)), '/')
    stopifnot(max(abs(lambdas-c(alpha^2/(6*pi^2),beta^2/(8*pi^2)))) < 1e-8,
              min(abs(diag(crossprod(phi, phi_exact)))) > .999999)
  }
  invgrid <- seq(0, 1, length.out=inverse_M)
  inverse <- matrix(vector('list', 15), 3, 5)
  min_density <- Inf
  for (g in 1:3) for (r in 1:5) {
    d <- 1 + as.vector(U(invgrid) %*% (theta[g,]*a[r,]))
    min_density <- min(min_density, min(d))
    if (min(d) < -1e-12) stop('Negative density in DGP')
    F <- invgrid + as.vector(H(invgrid) %*% (theta[g,]*a[r,]))
    F[1] <- 0; F[length(F)] <- 1
    if (any(diff(F) < -1e-12)) stop('Non-monotone CDF')
    inverse[[g,r]] <- approxfun(cummax(F), invgrid, ties='ordered', rule=2)
  }
  structure(list(grid=grid,w=w,sw=rep(sqrt(w),times=5L),
                 p=p,means=means,phi=phi,lambda=lambdas,A=A,
                 inverse=inverse, min_density=min_density,
                 alpha=alpha,beta=beta), class='smpm_config')
}

smpm_simulate <- function(cfg, n, depth=1000L, store_profiles=FALSE) {
  stopifnot(n>=3, length(depth)==1L, depth>=1, depth==as.integer(depth))
  y <- sample.int(3L,n,replace=TRUE,prob=cfg$p)
  counts <- tabulate(y, nbins=3L)
  if (any(counts==0)) stop('Empty phenotype group')
  d <- 5L*length(cfg$grid)
  sums <- matrix(0,3L,d)
  X <- if (store_profiles) matrix(0,n,d) else NULL
  for (i in seq_len(n)) {
    row <- numeric(d)
    for (r in 1:5) {
      events <- sort(cfg$inverse[[y[i],r]](runif(depth)))
      idx <- ((r-1L)*length(cfg$grid)+1L):(r*length(cfg$grid))
      row[idx] <- findInterval(cfg$grid, events, rightmost.closed=TRUE)/depth * cfg$sw[idx]
    }
    sums[y[i],] <- sums[y[i],] + row
    if (store_profiles) X[i,] <- row
  }
  group_means <- sweep(sums,1,counts,'/')
  list(y=y, counts=counts, group_means=group_means, X=X)
}

smpm_fit <- function(sim, reference=NULL) {
  p <- sim$counts/sum(sim$counts)
  m <- sim$group_means
  z <- colSums(m*p)
  V <- sweep(m,2,z,'-')
  A <- sweep(V,1,sqrt(p),'*')
  ee <- eigen(tcrossprod(A), symmetric=TRUE)
  lam <- pmax(ee$values[1:2],0)
  if (any(lam < 1e-14)) stop('Near-zero supervised eigenvalue')
  phi <- t(A)%*%ee$vectors[,1:2,drop=FALSE]
  phi <- sweep(phi,2,sqrt(lam),'/')
  out <- list(lambda=lam,phi=phi,group_means=m,group_probs=p)
  if (!is.null(reference)) {
    alignment <- abs(diag(crossprod(phi,reference$phi)))
    alignment <- pmin(1,pmax(0,alignment))
    out$angle <- acos(alignment)
    out$phi_error <- sqrt(pmax(0,2-2*alignment))
    # ||A'A - B'B||_F^2 = ||AA'||_F^2 + ||BB'||_F^2 - 2||AB'||_F^2
    cross <- A%*%t(reference$A)
    hs2 <- sum(tcrossprod(A)^2)+sum(tcrossprod(reference$A)^2)-2*sum(cross^2)
    out$hs_error <- sqrt(max(0,hs2))
  }
  if (!is.null(sim$X)) {
    idx <- sim$y
    a <- V%*%phi
    ar <- a[idx,,drop=FALSE]
    residual <- sim$X - m[idx,,drop=FALSE]
    psi <- ar^2+2*ar*(residual%*%phi)
    psi <- sweep(psi,2,lam,'-')
    out$influence <- psi
    out$sigma2 <- colMeans(sweep(psi,2,colMeans(psi),'-')^2)
  }
  out
}
