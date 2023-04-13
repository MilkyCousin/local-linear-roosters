source("../core/estimator.R")

cv.integrated.least.squares <- function(x, y, kernel, h, W.left, W.right, m, is.loc.lin=TRUE)
{
  n <- length(x)
  n.half <- floor(n/2)
  
  a.half.left <- W.left[,m]
  a.half.right <- W.right[,m]
  
  #a.half.left <- cormid(W.left[,m], order(x[1:n.half]))
  #a.half.right <- cormid(W.right[,m], order(x[(n.half+1):n]))

  # Підрахунок оцінки для g.m на основі всіх спостережень
  g.m.full <- sapply((n.half+1):n, function(j) {
    loc.lin.est(x[j], x[(n.half+1):n], y[(n.half+1):n], kernel, h, a.half.right, is.loc.lin)
  })
  # Підрахунок оцінки для g.m на основі всіх спостережень, окрім j-го
  g.m.out <- sapply(1:n.half, function(j) {
    loc.lin.est(x[j], x[1:n.half][-j], y[1:n.half][-j], kernel, h, a.half.left[-j], is.loc.lin)
  })
  
  #I1 <- sum(a.half.left * g.m.out * g.m.out)
  #I2 <- sum(a.half.left * g.m.out * g.m.full)
  #sum.components <- I1 - 2 * I2
  
  sum.components <- a.half.left * g.m.out * (g.m.out - 2 * g.m.full)
  result <- sum(sum.components)
  result
}

cv.least.squares <- function(x, y, kernel, h, a.m, is.loc.lin=TRUE)
{
  n <- length(x)
  jdx <- 1:n
  # Підрахунок оцінки для g.m на основі всіх спостережень, окрім j-го
  g.m.out <- sapply(jdx, function(j) {
    loc.lin.est(x[j], x[-j], y[-j], kernel, h, a.m[-j], is.loc.lin)
  })
  sum(a.m * (y - g.m.out)^2)
}