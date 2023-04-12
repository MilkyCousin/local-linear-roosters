## Завантаження попередніх розробок
source("../core/mixvconc.R")
source("../core/estimator.R")
source("../core/functions.R")

# Функція для наближеного підрахунку таких характеристик якості збіжності як
# 1. Зсунутість оцінки: Bias(\hat{\theta}) = E[\hat{\theta} - \theta]
# 2. Дисперсія оцінки: D[\hat{\theta}]
# Підрахунок оцінки робиться у заданій точці x0.
# Вхід:
# n.samples - кількість повторних вибірок з заданого розподілу
# x0 - точка, в якій оцінюється значення невідомих функцій зв'язку
# g.real - фабрика справжніх функцій зв'язку для кожної компоненти суміші
# params.d - параметри для генерування вибірки з суміші з залежністю:
# $x, $conc.matr, $x.distr, $e.distr, $g.func відповідають параметрам
# params.e - параметри для побудови непараметричної оцінки функції зв'язку:
# $kern - ядро (функція), 
# $band - параметр згладжування (функція від n)
# функції gen.sampl.mixt.regr(...)
# is.ll -- TRUE, якщо рахувати лок-лін оцінку, інакше Надарая-Ватсона
# seed.val - зернина для імітаційного моделювання
# Вихід:
# calculated.measures -- Таблиця з підрахованими характеристиками
imit.modeling.scheme.pointwise <- function(
    n.samples, x0, g.real, params.d, params.e, is.ll=TRUE, seed.val=0
)
{
  # Фіксація зернини для генератора псевдовипадкових чисел
  set.seed(seed.val)
  
  # Побудова функції для створення вибірки з моделі суміші
  gen.distr.data.current <- function()
  {
    to.return <- gen.sample.mixt.regr(
      params.d$n, 
      params.d$conc.matr, 
      params.d$x.distr, 
      params.d$e.distr, 
      params.d$g.func
    )
    to.return
  }
  
  # Кількість компонент в моделі суміші
  M <- ncol(params.d$conc.matr)
  
  # Обчислення мінімаксних вагових коефіцієнтів
  W.matr <- minmax.weights(params.d$conc.matr)
  
  # Побудова функції для підрахунку оцінки невідомої функції зв'язку
  stat.func.current <- function(m)
  {
    factory.func <- function(src.sample)
    {
      estimated.value <- loc.lin.est(
        x0, 
        src.sample$x, src.sample$y, 
        params.e$kern,
        params.e$band(params.d$n), 
        W.matr[,m],
        is.ll
      )
      estimated.value
    }
    factory.func
  }
  
  # Для кожної функції зв'язку, підраховуємо відповідні значення оцінок
  estimated.values <- calculate.multiple.statistics(
    gen.distr.data=gen.distr.data.current, 
    n.copies=n.samples, 
    stat.func=stat.func.current, 
    n.stats=M
  )
  
  # Вектор зі справжніх значень невідомих функцій зв'язку в точці
  theta.m <- sapply(1:M, function(m) { g.real(m)(x0) })
  
  # Наближений підрахунок зсунутості, дисперсії та MSE оцінок
  bias.m <- apply(estimated.values - theta.m, 1, mean)
  var.m <- apply(estimated.values, 1, var)
  # mse.m <- var.m + bias.m^2

  # Таблиця обчислених характеристик, яку має повернути функція
  calculated.measures <- data.frame(
    bias.val = bias.m,
    var.val = var.m
  )
  
  calculated.measures
}


# Функція для наближеного підрахунку такої характеристики якості збіжності як
# 1. MISE оцінки: E[\int (\hat{\theta}(x) - \theta(x))^2 dx]
# Підрахунок оцінки робиться у заданій точці x0.
# Вхід:
# n.samples - кількість повторних вибірок з заданого розподілу
# integr.param - параметри для підрахунку ISE
# $lower - нижня межа інтегрування, $upper - верхня межа інтегрування
# $subdiv - кількість поділів в інтегруванні
# $err.ignore - чи ігнорувати помилки в integrate(...)
# g.real - фабрика справжніх функцій зв'язку для кожної компоненти суміші
# params.d - параметри для генерування вибірки з суміші з залежністю:
# $x, $conc.matr, $x.distr, $e.distr, $g.func відповідають параметрам
# params.e - параметри для побудови непараметричної оцінки функції зв'язку:
# $kern - ядро (функція), 
# $band - параметр згладжування (функція від n)
# функції gen.sampl.mixt.regr(...)
# is.ll -- TRUE, якщо рахувати лок-лін оцінку, інакше Надарая-Ватсона
# seed.val - зернина для імітаційного моделювання
# Вихід:
# calculated.measures -- Таблиця з підрахованими характеристиками
imit.modeling.scheme.integration <- function(
    n.samples, integr.param, g.real, params.d, params.e, is.ll=TRUE, seed.val=0
)
{
  # Фіксація зернини для генератора псевдовипадкових чисел
  set.seed(seed.val)
  
  # Побудова функції для створення вибірки з моделі суміші
  gen.distr.data.current <- function()
  {
    to.return <- gen.sample.mixt.regr(
      params.d$n, 
      params.d$conc.matr, 
      params.d$x.distr, 
      params.d$e.distr, 
      params.d$g.func
    )
    to.return
  }
  
  # Кількість компонент в моделі суміші
  M <- ncol(params.d$conc.matr)
  
  # Обчислення мінімаксних вагових коефіцієнтів
  W.matr <- minmax.weights(params.d$conc.matr)
  
  # Побудова функції для підрахунку MISE оцінки невідомої функції зв'язку
  stat.func.current <- function(m)
  { 
    factory.func <- function(src.sample)
    {
      # Оцінка функції зв'язку
      theta.estim.univar <- function(t)
      {
        to.return <- loc.lin.est(
          t, 
          src.sample$x, src.sample$y, 
          params.e$kern, params.e$band(params.d$n), 
          W.matr[,m],
          is.ll
        )
        to.return
      }
      
      theta.estim <- function(t) { sapply(t, theta.estim.univar) }
      
      # Справжня функція зв'язку
      theta.real <- function(t)
      {
        to.return <- g.real(m)(t)
        to.return
      }
      
      # Підрахунок ISE
      estimated.value <- integrate(
        f = function(t) { (theta.estim(t) - theta.real(t))^2 },
        lower = integr.param$lower,
        upper = integr.param$upper,
        subdivisions = integr.param$subdiv,
        stop.on.error = integr.param$err.ignore
      )$value
      
      # Повертаємо ISE
      estimated.value
    }
    factory.func
  }
  
  # Для кожної функції зв'язку, підраховуємо ISE на повторних вибірках
  estimated.values <- calculate.multiple.statistics(
    gen.distr.data=gen.distr.data.current, 
    n.copies=n.samples, 
    stat.func=stat.func.current, 
    n.stats=M
  )
  
  # Наближений підрахунок MISE
  mise.m <- apply(estimated.values, 1, mean)
  
  # Таблиця обчислених характеристик, яку має повернути функція
  calculated.measures <- data.frame(mise.val = mise.m)
  
  calculated.measures
}