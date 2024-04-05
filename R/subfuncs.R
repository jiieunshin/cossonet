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

  if (type == "poly") {
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

  if (type == "spline") {
    K = 0
    for (d in 1:p) {
      K_temp = spline_kernel(x[, d, drop = FALSE], y[, d, drop = FALSE])
      K = K + K_temp$K1 + K_temp$K2
    }
  }

  if (type == "linear") {
    K = tcrossprod(x, y)
  }

  if (type == "anova_gaussian") {
    K = 0
    for (d in 1:p) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, type = "gaussian", kparam = kparam)
      K = K + K_temp
    }
  }

  return(K)
}

make_anovaKernel = function(x, y, type, kparam)
{
  x = as.matrix(x)
  y = as.matrix(y)
  dimx = ncol(x)

  # calculate anova kernels for main effects
  if (type == "spline") {
    # assign the number of anova kernels
    numK = 2 * dimx
    # list of kernel matrices
    anova_kernel = vector(mode = "list", numK)
    # list of kernel coordinate indices
    kernelCoord = vector(mode = "list", numK)
    index = 0

    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep="")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep="")
    }

  } else if (type == 'spline2') {
    numK = (2 * dimx) + (2 * dimx * (2 * dimx - 1) / 2 - dimx)
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    # main effects
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep = "")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep = "")
    }
    # two-way interactions
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A_linear = as.matrix(anova_kernel[[2 * i - 1]])
        A_smooth = as.matrix(anova_kernel[[2 * i]])
        B_linear = as.matrix(anova_kernel[[2 * j - 1]])
        B_smooth = as.matrix(anova_kernel[[2 * j]])
        anova_kernel[[index]] = A_linear * B_linear
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_linear * B_smooth
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " smooth", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_linear
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_smooth
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " smooth", sep = "")
      }
    }
  } else if (type == "spline-t") {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
  } else if (type == 'spline-t2') {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0

    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
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
      anova_kernel[[index]] = kernelMatrix(A, B, type, kparam)
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
  } else {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    for (d in 1:dimx) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      anova_kernel[[d]] = kernelMatrix(A, B, type, kparam)
      kernelCoord[[d]] = paste("x", d, sep = "")
    }
  }
  return(list(x = x, K = anova_kernel, coord = kernelCoord, numK = numK, kernel = type, kparam = kparam))
}

# is used to combine anova kernel matrices with weights determined by theta values. The default theta vector is the vector of ones.
combine_kernel = function(anova_kernel, theta = rep(1, anova_kernel$numK))
{
  K = 0
  for (d in 1:anova_kernel$numK) {
    K = (K + theta[d] * anova_kernel$K[[d]])
  }
  return(K)
}

rescale = function (x)
{
  if (length(unique(x)) > 6)
    return((x - min(x))/(max(x) - min(x)))
  else return(x)
}

wsGram = function (Gramat, mscale)
{
  n1 <- dim(Gramat)[1]
  n2 <- dim(Gramat)[2]
  d <- dim(Gramat)[3]
  KK <- matrix(0, n1, n2)
  for (j in 1:d) KK = KK + mscale[j] * Gramat[, , j]
  return(KK)
}

rescale = function (x)
{
  if (length(unique(x)) > 6)
    return((x - min(x))/(max(x) - min(x)))
  else return(x)
}

cvsplitID = function (n, folds)
{
  fsize <- floor(n/folds)
  splits <- fsize * rep(1, folds)
  nextra <- n - folds * fsize
  if (nextra > 0) {
    splits[1:nextra] <- splits[1:nextra] + 1
  }
  randid <- sample(1:n, n)
  IDmat <- matrix(NA, ncol = folds, nrow = ceiling(n/folds))
  IDmat[, 1] <- randid[1:splits[1]]
  for (i in 2:folds) {
    tempid <- randid[(cumsum(splits)[i - 1] + 1):(cumsum(splits)[i])]
    length(tempid) <- ceiling(n/folds)
    IDmat[, i] <- tempid
  }
  return(IDmat)
}

rescale_theta = function(theta, scale = FALSE){
  d = length(theta)
  if(sum(theta == 0) == d){
    theta = rep(1e-10, d)
  } else{
    if(scale) theta = theta / sd(theta)
    if(!scale) theta = theta
  }
  return(theta)
}

KLD = function(y, fhat, obj){
  if(obj$family == "gaussian") B = function(x) x
  if(obj$family == "binomial") B = function(x) log(exp(x) + 1)
  if(obj$family == "poisson") B = function(x) exp(x)

  return(- y * fhat + B(fhat))
}
