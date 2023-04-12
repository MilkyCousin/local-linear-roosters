  ## Завантаження попередніх розробок
  source("./imitation.R")
  source("./crossvalidation.R")
  
  ### Застосування розробок на імітаційне моделювання на конкретних прикладах
  
  ## Параметри оцінки
  # Вибір параметрів непараметричної оцінки
  k.epan <- function(t) { 0.75 * (abs(t) < 1) * (1 - t^2) }
  h.func <- function(k) { k^(-1/5) }
  params.e.current <- list(kern = k.epan, band = h.func)
  
  ## Параметри розподілу
  # Кількість повторних вибірок
  given.n.samples <- 10000
  
  # Обсяг вибірки: 100, 500, 1000, 2500, 5000, 10000
  given.n <- c(100, 500, 1000, 2500, 5000, 10000)
  
  # Який саме варіант функцій зв'язку та похибок розглядається?
  
  # 1 - дві параболи
  # 2 - розрив
  # 3 - тригонометрія
  # 4 - близькі з перетином
  # 5 - модуль та парабола
  example.number.func <- 1
  g.func.selected <- g.func.m(example.num = example.number.func)
  
  # 1 - N(0, 0.0025), 2 - N(0, 1.25), 3 - T(5)
  example.number.errors <- 1
  e.distr.m.selected <- e.distr.m(example.num = example.number.errors)
  
  ## Реалізація крос-валідації
  n.current <- 10000
  
  # Матриця концентрацій
  given.conc.matr <- gen.two.comp.prob(n.current)
  
  # Вибір параметрів розподілу в моделі суміші
  #params.d.current <- list(
  #  n = n.current, 
  #  conc.matr = given.conc.matr, 
  #  x.distr = x.distr.m,
  #  e.distr = e.distr.m.selected,
  #  g.func = g.func.selected
  #)
  
  gen.sample <- gen.sample.mixt.regr(
    n=n.current, 
    conc.matr=given.conc.matr, 
    x.distr=x.distr.m, 
    e.distr=e.distr.m.selected, 
    g.func=g.func.selected
  )
  
  # Обчислення мінімаксних вагових коефіцієнтів
  W.matr <- minmax.weights(given.conc.matr)
  
  # Параметри на крос-валідацію
  h.start <- h.func(n.current)
  h.values <- seq(0.001 * h.start, 10 * h.start, 0.025)
  
  W.matr.1 <- cormid(W.matr[,2], order(gen.sample$x))
  cv.1 <- function(u) 
  {
    #cv.functional.1(gen.sample$x, gen.sample$y, k.epan, u, W.matr.1)
    cv.functional.2(gen.sample$x, gen.sample$y, k.epan, u, given.conc.matr, 2)
  }
  
  h0 <- n.current^(-1/5)
  h.cv.1 <- optim(h0, cv.1, method="Brent",
                  lower=0.001 * h0, upper=3 * h0)$par
  #h.values[which.min(cv.1.values)]
  
  print("diagram")
  cv.1.values <- sapply(h.values, cv.1)
  
  plot(h.values, cv.1.values, type='l'); grid()
  abline(v=h.cv.1, col='red')
  abline(v=h.start, col='blue')
  
  # Діаграма розсіювання даних
  pal <- c('red', 'blue')
  plot(gen.sample$x, gen.sample$y, col=pal[gen.sample$m],
       xlab='x', ylab='y', cex=0.5,
       main=paste(
         "n = ", n.current, sep=""
       ))
  grid()
  
  # Побудова оцінок функцій зв'язку
  stat.func.current <- function(m)
  { 
    factory.func <- function(t)
    {
      estimated.value <- loc.lin.est(
        t, gen.sample$x, gen.sample$y, k.epan, h.cv.1, W.xy[,m], TRUE
      )
      estimated.value
    }
    
    factory.func.vec <- function(t) { sapply(t, factory.func) }
    factory.func.vec
  }
  
  # Побудова графіків справжніх функцій зв'язку
  line.thickness <- 3
  
  curve(g.func.selected(2)(x), col="green", add=T, lwd=line.thickness, lty=2)
