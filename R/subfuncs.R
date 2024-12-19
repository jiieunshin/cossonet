make_kernel = function (x, y, type)
{
  n1 <- nrow(x)
  n2 <- nrow(y)
  d <- ncol(x)
  K <- array(0, c(n1, n2, d))
  for (j in 1:d) {
    K[, , j] <- kernelMatrix(x, y, type)
  }
  return(K)
}

spline_kernel = function(x, y)
{
  x = as.matrix(x)
  y = as.matrix(y)
  K1x = (x - 1 / 2)
  K1y = (y - 1 / 2)
  K2x = (K1x^2 - 1 / 12) / 2
  K2y = (K1y^2 - 1 / 12) / 2
  ax = x %x% matrix(1, 1, nrow(y))
  ay = y %x% matrix(1, 1, nrow(x))
  b = abs(ax - t(ay))
  K1 = K1x %x% t(K1y)
  K2 = K2x %x% t(K2y) - ((b - 1 / 2)^4 - (b - 1 / 2)^2 / 2 + 7 / 240) / 24
  list(K1 = K1, K2 = K2)
}

cat_kernel = function(x, y)
{
  x = as.matrix(x)
  y = as.matrix(y)

  n1 <- length(x)
  n2 <- length(y)
  x <- rep(x, times = n2)
  y <- rep(y, each = n1)
  L <- length(unique(c(x, y)))
  K <- matrix(L * (x == y) - 1, n1, n2)
  return(K)
}

kernelMatrix = function(x, y, type, kparam = 1.0) {

  x = as.matrix(x)
  y = as.matrix(y)
  p = ncol(x)

  if (ncol(x) == 0) {
    x = matrix(0, nrow = nrow(x), ncol = 1)
  }

  if (ncol(y) == 0) {
    y = matrix(0, nrow = nrow(y), ncol = 1)
  }

  if (type == "poly" | type == "poly2") {
    K = (x %*% t(y) + 1.0)^kparam
  }

  if(type == "gaussian" | type == "gaussian2") {
    normx = rowSums(x^2)
    normy = rowSums(y^2)
    temp = x %*% t(y)
    temp = (-2.0 * temp) + outer(normx, rep(1.0, nrow(y)), "*") + outer(rep(1.0, nrow(x)), normy, "*")
    K = exp(-temp * kparam)
    # obj = kernelMatrix(rbfdot(sigma = kparam), x, y)
  }

  if (type == "spline" | type == "spline2") {
    K = 0
    for (d in 1:p) {
      K_temp = spline_kernel(x[, d, drop = FALSE], y[, d, drop = FALSE])
      K = K + K_temp$K1 + K_temp$K2
    }
  }

  if (type == "linear" | type == "linear2") {
    K = tcrossprod(x, y)
  }

  return(K)
}

make_anovaKernel = function(x, y, type, kparam, scale)
{

  # if (length(unique(c(A, B))) <= 6)
  #   K_temp <- cat_kernel(A, B)
  # else K_temp <- spline_kernel(A, B)

  x = as.matrix(x)
  y = as.matrix(y)
  dimx = ncol(x)

  # calculate anova kernels for two-way interactions
  if (type == "spline") {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      if (length(unique(c(A, B))) <= 6){
        K_temp <- cat_kernel(A, B)
        anova_kernel[[index]] = K_temp
      } else{
        K_temp = spline_kernel(A, B)
        anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      }
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
  } else if (type == 'spline2') {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0

    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      if (length(unique(c(A, B))) <= 6){
        K_temp <- cat_kernel(A, B)
        anova_kernel[[index]] = K_temp
      } else{
        K_temp = spline_kernel(A, B)
        anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      }
      kernelCoord[[index]] = paste("x", d, sep = "")
    }

    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  } else if (type == "gaussian2") {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0

    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      if (length(unique(c(A, B))) <= 6){
        K_temp <- cat_kernel(A, B)
        anova_kernel[[index]] = K_temp
      } else{
        anova_kernel[[index]] = kernelMatrix(A, B, type, kparam)
      }
      kernelCoord[[index]] = paste("x", d, sep = "")
    }

    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  } else { # calculate anova kernels for main effects
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    for (d in 1:dimx) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      if (length(unique(c(A, B))) <= 6){
        K_temp <- cat_kernel(A, B)
        anova_kernel[[d]] = K_temp
      } else{
        anova_kernel[[d]] = kernelMatrix(A, B, type, kparam)
      }
      kernelCoord[[d]] = paste("x", d, sep = "")
    }
  }


  return(list(x = x, K = anova_kernel, coord = kernelCoord, numK = numK, kernel = type, kparam = kparam))
}

rescale = function (x)
{
  if (length(unique(x)) > 6)
    return((x - min(x))/(max(x) - min(x)))
  else return(x)
}

combine_kernel = function (Gramat, mscale)
{
  n1 <- dim(Gramat)[1]
  n2 <- dim(Gramat)[2]
  d <- dim(Gramat)[3]
  KK <- matrix(0, n1, n2)
  for (j in 1:d) KK = KK + mscale[j] * Gramat[, , j]
  return(KK)
}

rescale_theta = function (x)
{
  d = length(x)
  if(sum(x == 0) == d) x = rep(1e-10, d)
  return(x)
}

cvsplitID = function (n, folds, y, family)
{
  fsize <- floor(n/folds)
  splits <- fsize * rep(1, folds)
  nextra <- n - folds * fsize
  if (nextra > 0) {
    splits[1:nextra] <- splits[1:nextra] + 1
  }

  if(family != "binomial"){
    randid <- sample(1:n, n)
    IDmat <- matrix(NA, ncol = folds, nrow = ceiling(n/folds))
    IDmat[, 1] <- randid[1:splits[1]]
    for (i in 2:folds) {
      tempid <- randid[(cumsum(splits)[i - 1] + 1):(cumsum(splits)[i])]
      length(tempid) <- ceiling(n/folds)
      IDmat[, i] <- tempid
    }
  }

  if(family == "binomial"){
    if(is.null(y)) stop("The input of y is essential.")

    # Separate indices for 0s and 1s
    idx_0 <- which(y == 0)
    idx_1 <- which(y == 1)

    n0 <- length(idx_0)
    n1 <- length(idx_1)

    # Compute fold sizes for each class
    fsize_0 <- floor(n0 / folds)
    fsize_1 <- floor(n1 / folds)

    splits_0 <- fsize_0 * rep(1, folds)
    splits_1 <- fsize_1 * rep(1, folds)

    nextra_0 <- n0 - folds * fsize_0
    nextra_1 <- n1 - folds * fsize_1

    if (nextra_0 > 0) splits_0[1:nextra_0] <- splits_0[1:nextra_0] + 1
    if (nextra_1 > 0) splits_1[1:nextra_1] <- splits_1[1:nextra_1] + 1

    randid_0 <- sample(idx_0, n0)
    randid_1 <- sample(idx_1, n1)

    IDmat <- matrix(NA, ncol = folds, nrow = ceiling(n / folds))

    # Assign 0s and 1s to folds
    for (i in 1:folds) {
      if(i == 1){
        tempid_0 <- randid_0[1:(cumsum(splits_0)[i])]
        tempid_1 <- randid_1[1:(cumsum(splits_1)[i])]
      } else{
        tempid_0 <- randid_0[(cumsum(splits_0)[i - 1] + 1):(cumsum(splits_0)[i])]
        tempid_1 <- randid_1[(cumsum(splits_1)[i - 1] + 1):(cumsum(splits_1)[i])]
      }

      tempid <- c(tempid_0, tempid_1)
      length(tempid) <- ceiling(n / folds)
      IDmat[, i] <- tempid
    }
  }
  return(IDmat)
}


KL = function(f, mu, obj){
  if(obj$family == "gaussian") B = f^2/2
  if(obj$family == "binomial") B = log(1 + exp(f))
  if(obj$family == "poisson") B = exp(f)

  return(mean(-(mu * f) + B))
}


# SKL = function(f, fhat){
#   return(mean(f * log(f / fhat)))
# }
