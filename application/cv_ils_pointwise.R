### Крос-валідація на основі середньквадратичної похибки. Поточкова якість

## Завантаження попередніх розробок
source("./imitation.R")
source("./crossvalidation.R")

# Вибір параметрів непараметричної оцінки
k.epan <- function(t) { 0.75 * (abs(t) < 1) * (1 - t^2) }
h.func <- function(k, ...) { k^(-1/5) }
k.func <- k.epan

# Кількість повторних вибірок
given.n.samples <- 1000

# Обсяг вибірки: 100, 500, 1000, 2500, 5000, 10000
#given.n <- c(2500) 
given.n <- c(100, 500, 1000)#, 2500, 5000, 10000)

# Точка, в якій оцінюються функції зв'язку
given.x0 <- 0.5

# Який саме варіант функцій зв'язку та похибок розглядається?

# 1 - дві параболи
# 2 - розрив
# 3 - тригонометрія
# 4 - близькі з перетином
# 5 - модуль та парабола
example.number.func <- 1

# Функції зв'язку
g.func.selected <- g.func.m(example.num = example.number.func)
# Справжні значення в точках
theta.1 <- g.func.selected(1)(given.x0)
theta.2 <- g.func.selected(2)(given.x0)

# 1 - N(0, 0.0025), 2 - N(0, 1.25), 3 - T(5)
example.number.errors <- 2
e.distr.m.selected <- e.distr.m(example.num = example.number.errors)

print(paste("Приклад функції зв'язку: ", example.number.func))
print(paste("Приклад розподілу похибок: ", example.number.errors))

## Імітаційне моделювання

# Для кожного обсягу вибірки проводимо імітаційний експеримент
for(n.current in given.n)
{
  pointwise.results.1 <- numeric(given.n.samples)
  pointwise.results.2 <- numeric(given.n.samples)
  
  print(paste("- Обробка для n =", n.current))
  
  # Створюємо матрицю концентрнацій
  given.conc.matr <- gen.two.comp.prob(n.current)
  
  # Створюємо мінімаксні вагові коефіцієнти
  n.half <- floor(n.current / 2)
  
  W.matr <- minmax.weights(given.conc.matr)
  W.left <- minmax.weights(given.conc.matr[1:n.half,])
  W.right <- minmax.weights(given.conc.matr[(n.half+1):n.current,])
  
  # Генеруємо B := given.n.samples копій вибірки із заданим розподілом
  for(b in 1:given.n.samples)
  {
    if(b %% 50 == 0)
    {
      print(paste("Пройшло", b, "Етапів"))
    }
    # Генеруємо вибірку з регресійної моделі суміші
    xym.sample <- gen.sample.mixt.regr(
      n.current, given.conc.matr, 
      x.distr.m, e.distr.m.selected, 
      g.func.selected
    )
    
    # Обираємо параметр згладжування на основі середньоквадратричної CV
    h.init <- n.current^(-1/5)
    
    # Функціонали CV, навантажені на перші дві компоненти (інших не буде)
    cv.f.1 <- function(u) 
    { 
      cv.integrated.least.squares(xym.sample$x, xym.sample$y, k.func, u, W.left, W.right, 1)
    }
    cv.f.2 <- function(u) 
    { 
      cv.integrated.least.squares(xym.sample$x, xym.sample$y, k.func, u, W.left, W.right, 2)
    }
    
    # Вибір параметрів згладжування оцінки
    h.cv.1 <- optimise(
      f=cv.f.1, lower=0.01*h.init, upper=5*h.init
    )$minimum
    h.cv.2 <- optimise(
      f=cv.f.2, lower=0.01*h.init, upper=5*h.init
    )$minimum
    
    # Підрахунок оцінок в заданих точках
    estim.1 <- loc.lin.est(
      given.x0, xym.sample$x, xym.sample$y, k.func, h.cv.1, W.matr[,1], TRUE
    )
    estim.2 <- loc.lin.est(
      given.x0, xym.sample$x, xym.sample$y, k.func, h.cv.2, W.matr[,2], TRUE
    )
    
    # Записуємо в масив
    pointwise.results.1[b] <- estim.1
    pointwise.results.2[b] <- estim.2
  }
  
  # Результати моделювання:
  print("Перша компонента:")
  print(c("Зміщення", "Дисперсія"))
  print(c(mean(pointwise.results.1 - theta.1), var(pointwise.results.1)))
  
  print("Друга компонента:")
  print(c("Зміщення", "Дисперсія"))
  print(c(mean(pointwise.results.2 - theta.2), var(pointwise.results.2)))
}