### Застосування розробок на імітаційне моделювання на конкретних прикладах

## Завантаження попередніх розробок
source("./imitation.R")
source("./crossvalidation.R")

# Вибір параметрів непараметричної оцінки
k.epan <- function(t) { 0.75 * (abs(t) < 1) * (1 - t^2) }
h.func <- function(k, ...) { k^(-1/5) }
params.e.current <- list(kern = k.epan, band = h.func)

# Кількість повторних вибірок
given.n.samples <- 1000

# Обсяг вибірки: 100, 500, 1000, 2500, 5000, 10000
given.n <- c(100, 500, 1000, 2500, 5000, 10000)

# Точки, в яких оцінюються функції зв'язку
given.x0 <- c(0, 0.5, 1)

# Який саме варіант функцій зв'язку та похибок розглядається?

# 1 - дві параболи
# 2 - розрив
# 3 - тригонометрія
# 4 - близькі з перетином
# 5 - модуль та парабола
example.number.func <- 1
g.func.selected <- g.func.m(example.num = example.number.func)

# 1 - N(0, 0.0025), 2 - N(0, 1.25), 3 - T(5)
example.number.errors <- 2
e.distr.m.selected <- e.distr.m(example.num = example.number.errors)

print(paste("Приклад функції зв'язку: ", example.number.func))
print(paste("Приклад розподілу похибок: ", example.number.errors))

## Імітаційне моделювання

pointwise.results <- list()

for(x0.current in given.x0)
{
  print(paste("Обробка для x0 =", x0.current))
  
  # Ініціалізація таблиць
  table.ll <- data.frame(
    n = numeric(),
    bias.val.1 = numeric(),
    var.val.1 = numeric(),
    bias.val.2 = numeric(),
    var.val.2 = numeric()
  )
  table.nw <- data.frame(
    n = numeric(),
    bias.val.1 = numeric(),
    var.val.1 = numeric(),
    bias.val.2 = numeric(),
    var.val.2 = numeric()
  )
  
  # Наповнення таблиць
  for(n.current in given.n)
  {
    print(paste("- Обробка для n =", n.current))
    
    # Матриця концентрацій
    given.conc.matr <- gen.two.comp.prob(n.current)
    
    # Вибір параметрів розподілу в моделі суміші
    params.d.current <- list(
      n = n.current, 
      conc.matr = given.conc.matr, 
      x.distr = x.distr.m,
      e.distr = e.distr.m.selected,
      g.func = g.func.selected
    )
    
    # Цей етап можна теж оптимізувати. Одне моделювання замість двох
    imit.model.res.ll <- imit.modeling.scheme.pointwise(
      n.samples = given.n.samples, 
      x0 = x0.current, 
      g.real = g.func.selected,
      params.d = params.d.current,
      params.e = params.e.current,
      is.ll = TRUE
    )
    imit.model.res.nw <- imit.modeling.scheme.pointwise(
      n.samples = given.n.samples, 
      x0 = x0.current, 
      g.real = g.func.selected,
      params.d = params.d.current,
      params.e = params.e.current,
      is.ll = FALSE
    )
    
    # Збереження таблиць
    for.rows.ll <- cbind(
      n.current, imit.model.res.ll[1,], imit.model.res.ll[2,]
    )
    for.rows.nw <- cbind(
      n.current, imit.model.res.nw[1,], imit.model.res.nw[2,]
    )
    
    table.ll[nrow(table.ll) + 1, ] <- for.rows.ll
    table.nw[nrow(table.nw) + 1, ] <- for.rows.nw
  }
  
  # Додавання таблиць
  to.append <- list(
    x0 = x0.current,
    res.ll = table.ll,
    res.nw = table.nw
  )
  pointwise.results <- append(pointwise.results, to.append)
}

# Записуємо результати в текстовий файл
file.out <- paste('./examples/example', example.number.func, '/pointwise',
                  '-err', example.number.errors, 
                  '-func', example.number.func, 
                  '.txt', sep='')
sink(file.out)
print(pointwise.results)
sink()