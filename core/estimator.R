# Побудова оцінки локально-лінійної регресії в точці
# Вхід:
# x0 -- точка, в якій оцінюється невідома функція зв'язку
# x -- вектор регресорів
# y -- вектор відгуків
# kernel -- ядро непараметричної оцінки
# h -- параметр згладжування непараметричної оцінки
# a.m -- вектор вагових коефіцієнтів, що відповідає m-ій компоненті суміші
# is.loc.lin -- TRUE, якщо рахувати лок-лін оцінку, інакше Надарая-Ватсона
# Вихід:
# b -- значення оцінки в заданій точці (число)
loc.lin.est.1 <- function(x0, x, y, kernel, h, a.m, is.loc.lin = TRUE)
{
  w <- kernel((x0 - x) / h)
  x.centered <- x - x0
  
  My <- sum(w * a.m * y)
  Mw <- sum(w * a.m)
  
  if(is.loc.lin)
  {
    Mx <- sum(w * a.m * x.centered)
    Mxx <- sum(w * a.m * x.centered^2)
    Mxy <- sum(w * a.m * x.centered * y)
    
    k <- (Mxy * Mw - Mx * My) / (Mxx * Mw - Mx * Mx)
    b <- My / Mw - k * Mx / Mw
  }
  else
  {
    b <- My / Mw
  }
  
  b
}

loc.lin.est <- function(x0, x, y, kernel, h, a.m, is.loc.lin = TRUE)
{
  w <- kernel((x0 - x) / h)
  x.centered <- x - x0
  
  s1 <- sum(a.m * w * x.centered)
  s2 <- sum(a.m * w * x.centered * x.centered)
  wj <- (s2 - s1 * x.centered) * a.m * w
  
  res <- sum(wj * y) / sum(wj)
  res
}