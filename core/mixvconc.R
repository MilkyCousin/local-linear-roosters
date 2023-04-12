# Генерування випадкового вектора, рівномірно розподіленого на симплексі
# Вхід:
# d -- розмірність простору, в якому мітситься симплекс
# Вихід:
# v / sum(v) -- випадкова точка на симплексі
gen.simplex <- function(d)
{
  v <- rexp(d)
  v / sum(v)
}

# Генерування матриці ймовірностей змішувань (матриці концентрацій)
# за допомогою рівномірно розподілених векторів на симплексі
# Вхід:
# n -- кількість об'єктів
# m -- кількість компонент
# Вихід:
# P -- матриця концентрацій
gen.mix.prob <- function(n, m)
{
  P <- matrix(t(replicate(n, gen.simplex(m))), nrow = n, ncol = m)
  colnames(P) <- paste("Comp.", 1:m, sep="")
  P
}

# Генерування матриці концентрацій для двокомпонентної суміші
gen.two.comp.prob <- function(n)
{
  P <- matrix(nrow = n, ncol = 2)
  P[,1] <- 1:n / n
  P[,2] <- 1 - P[,1]
  P
}

# Обчислення мінімаксних вагових коефіцієнтів
# Вхід:
# p -- матриця концентрацій
# Вихід:
# p%*%gram.inv -- матриця мінімаксних вагових коефіцієнтів
minmax.weights <- function(p)
{
  p <- as.matrix(p)
  gram <- t(p)%*%p
  gram.inv <- solve(gram)
  p%*%gram.inv
}

# sup-корекція вагових коефіцієнтів
# Вхід:
# a -- початкові вагові коефіцієнти, 
# idx.sort -- антиранги,
# s -- наперед обчислені накопичені вагові коефіцієнти (якщо вказані)
# Вихід:
# b -- кориговані вагові коефіцієнти a
corsup <- function(a, idx.sort, s = NULL)
{
  if(is.null(s))
  {
    s <- c(0,cumsum(a[idx.sort]))
  }
  n <- length(a)
  s.plus <- pmin(1,cummax(s))
  b <- numeric(n)
  b[idx.sort] <- diff(s.plus)
  b
}

# inf-корекція вагових коефіцієнтів
# Вхід:
# a -- початкові вагові коефіцієнти, 
# idx.sort -- антиранги,
# s -- наперед обчислені накопичені вагові коефіцієнти (якщо вказані)
# Вихід:
# b -- кориговані вагові коефіцієнти a
corinf <- function(a, idx.sort, s = NULL)
{
  if(is.null(s))
  {
    s <- c(0,cumsum(a[idx.sort]))
  }
  n <- length(a)
  s.minus <- pmax(0,cummin(s[(n+1):1]))[(n+1):1]
  b <- numeric(n)
  b[idx.sort] <- diff(s.minus)
  b
}

# Комбінування sup- та inf- методів виправлення
# Вхід:
# a -- початкові вагові коефіцієнти, 
# idx.sort -- антиранги,
# p -- вплив sup-корекції, (1 - p) -- вплив inf-корекції
# Вихід:
# b -- кориговані вагові коефіцієнти a
cormid <- function(a, idx.sort, p = 0.5)
{
  n <- length(a)
  s <- c(0,cumsum(a[idx.sort]))
  p * corinf(a, idx.sort, s) + (1 - p) * corsup(a, idx.sort, s)
}

# Функція, яка генерує вибірку з моделі суміші зі змінним 
# концентраціями з залежністю
# Вхід:
# n -- кількість об'єктів
# conc.matr -- матриця концентрацій
# x.distr -- фабрика функцій, що реалізують умовні розподіли регресора
# e.distr -- фабрика функцій, що реалізують умовні розподіли похибки
# g.func -- фабрика функцій зв'язку кожної з компонент
# Вихід:
# to.return = list(x = numeric(n), y = numeric(n))
# -- змодельована вибірка з регресорами та відгуками обсягу n
gen.sample.mixt.regr <- function(n, conc.matr, x.distr, e.distr, g.func)
{
  # Кількість компонент суміші у моделі
  M <- ncol(conc.matr)
  
  # Генерування номерів компонент для n об'єктів вибірки
  idx <- sapply(
    1:n,
    function(j) {
      sample(1:M, size=1, replace=F, prob=conc.matr[j,])
    }
  )
  
  # Вектори регресорів та відгуків відповідно
  x.sample <- numeric(n)
  y.sample <- numeric(n)
  
  # На основі згенерованих номерів, генеруємо регресори
  # та похибки з заданими умовними розподілами
  for(m in 1:M)
  {
    n.m <- sum(idx == m)
    x.comp.m <- x.distr(m)(n.m)
    e.comp.m <- e.distr(m)(n.m)
    x.sample[idx == m] <- x.comp.m
    y.sample[idx == m] <- g.func(m)(x.comp.m) + e.comp.m
  }
  # Поєднуємо регресори з відгуками в один об'єкт, 
  # повертаємо об'єднане на виході
  to.return <- list(x=x.sample, y=y.sample, m=idx)
  to.return
}

# Моделює задану кількість повторних вибірок з обраного розподілу,
# на кожній з вибірок обчислюється статистика. В результі функція повертає
# масив зі значень статистики, обчислених за повторними вибірками
# Вхід:
# gen.distr.data -- конфігурація цільового розподілу даних, функція
# n.copies -- кількість повторних вибірок для створення
# stat.func -- статистики, які треба обчислити (фабрика функцій)
# n.stats -- кількість статистик для обчислення (краще поєднати з попереднім?)
# Вихід:
# stats.sample -- матриця зі значень обраних статистик (рядок -- оцінка)
calculate.multiple.statistics <- function(
    gen.distr.data, n.copies, stat.func, n.stats
)
{
  # Ініціалізація матриці для значень обраних оцінок
  stats.sample <- matrix(n.copies * n.stats, nrow = n.stats, ncol = n.copies)
  
  for(b in 1:n.copies)
  {
    if(b%%10==0)
    {
      print(paste("step:", b))
    }
    # Створюємо повторну вибірку
    current.sample <- gen.distr.data()
    
    for(m in 1:n.stats)
    {
      # Обчислюємо m-ту статистику за вибіркою
      stats.sample[m, b] <- stat.func(m)(current.sample)
    }
  }
  
  # Повертаємо масив зі значень статистик
  stats.sample
}