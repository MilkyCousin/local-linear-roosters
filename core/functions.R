### Приклади функцій зв'язків та розподілів регресорів, похибок. Можливо щось ще

## Програмна реалізація прикладів

# Розподіл регресора
x.distr.m <- function(m)
{
  to.return <- function(k) { runif(k, min=0, max=1) }
  to.return
}

# Розподіл похибки
#e.distr.m <- function(m)
#{
#  to.return <- function(k) { rnorm(k, mean=0, sd=0.025) }
#  to.return
#}

e.distr.m <- function(example.num)
{
  if(example.num == 1)
  {
    e.func <- function(m)
    {
      to.return <- function(k) { rnorm(k, mean=0, sd=sqrt(0.025)) }
      to.return
    }
  }
  else if(example.num == 2)
  {
    e.func <- function(m)
    {
      to.return <- function(k) { rnorm(k, mean=0, sd=sqrt(1.25)) }
      to.return
    }
  }
  else if(example.num == 3)
  {
    e.func <- function(m)
    {
      to.return <- function(k) { rt(k, df=5) }
      to.return
    }
  }
}

# Функції зв'язку
g.func.m <- function(example.num)
{
  if(example.num == 1)
  {
    g.func <- function(m)
    {
      to.return <- function(t) { (-1)^m * t * (1 - t) }
      to.return
    }
  }
  else if(example.num == 2)
  {
    g.func <- function(m)
    {
      if(m == 1)
      {
        to.return <- function(t) { sign(t - 0.7) * 3 * t * (t - 1) }
      }
      else
      {
        to.return <- function(t) { t + 1 }
      }
      to.return
    }
  }
  else if(example.num == 3)
  {
    g.func <- function(m)
    {
      to.return <- function(t) { sin(2 * pi * t + (m - 1) * pi / 2) }
      to.return
    }
  }
  else if(example.num == 4)
  {
    g.func <- function(m)
    {
      to.return <- function(t) { t + (m - 1) * sin(2 * pi * t) / 64 }
      to.return
    }
  }
  else if(example.num == 5)
  {
    g.func <- function(m)
    {
      to.return <- function(t) { abs(t - 0.5)^m }
      to.return
    }
  }
  g.func
}

## Інше

# Обчислення середньоквадратичної похибки (MSE) оцінки за вибіркою її значень
# Вхід:
# theta.estim -- масив значень оцінки для невідомого параметра theta
# theta.val -- невідомий параметр
# Вихід:
# to.return -- значення MSE
mse <- function(theta.estim, theta.val)
{
  to.return <- mean((theta.estim - theta.val)^2)
  to.return
}