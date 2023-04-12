## Завантаження попередніх розробок
source("../core/mixvconc.R")
source("../core/estimator.R")
source("../core/functions.R")

# Який саме приклад?
example.num.func <- 1
example.num.err <- 1
is.loc.lin.now <- FALSE

# Зернина
set.seed(1981 + example.num.func + example.num.err)

## Моделювання вибірки з двох компонент
n.given <- c(100, 500, 1000, 2500, 5000, 10000)

# Розмітка рисунку

png(filename = paste("../../report/pictures/", ifelse(is.loc.lin.now, "LL", "NW"),
                     "/func-", example.num.func, 
                     "-err-", example.num.err,
                     ".png",
                     sep = ""),
    width = 2400,
    height = 1800,
    res=256)

layout(
  mat=matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=TRUE),
  heights=c(1, 1, 0.1)
)

par(oma=c(0,0,3,0))
old.mai <- par("mai")

for(n in n.given)
{
  
  # Матриця концентрацій
  P <- gen.two.comp.prob(n)
  
  # Створення підвибірки
  
  g <- g.func.m(example.num = example.num.func)
  e.distr.m.selected <- e.distr.m(example.num = example.num.err)
  
  xy.dat <- gen.sample.mixt.regr(
    n=n, 
    conc.matr=P, 
    x.distr=x.distr.m, 
    e.distr=e.distr.m.selected, 
    g.func=g
  )
  
  # Діаграма розсіювання даних
  pal <- c('red', 'blue')
  plot(xy.dat$x, xy.dat$y, col=pal[xy.dat$m],
       xlab='x', ylab='y', cex=0.5,
       main=paste(
         "n = ", n, sep=""
       ))
  grid()
  
  # Побудова оцінок функцій зв'язку
  kern.xy <- function(t) { (abs(t) < 1) * (1 - t^2) }
  h.xy <- n^(-1/5)
  W.xy <- minmax.weights(P)
  
  stat.func.current <- function(m)
  { 
    factory.func <- function(t)
    {
      estimated.value <- loc.lin.est(
        t, xy.dat$x, xy.dat$y, kern.xy, h.xy, W.xy[,m], is.loc.lin.now
      )
      estimated.value
    }
    
    factory.func.vec <- function(t) { sapply(t, factory.func) }
    factory.func.vec
  }
  
  # Побудова графіків справжніх функцій зв'язку
  line.thickness <- 3
  
  curve(g(1)(x), col="green", add=T, lwd=line.thickness, lty=2)
  curve(g(2)(x), col="orange", add=T, lwd=line.thickness, lty=2)
  
  # Побудова графіків оцінок функцій зв'язку
  curve(stat.func.current(1)(x), col="darkgreen", add=T, lwd=line.thickness)
  curve(stat.func.current(2)(x), col="yellow", add=T, lwd=line.thickness)
}

# Створення легенди

par(mai=c(0,0,0,0))
plot.new()
legend("center", ncol=4,
       legend = c(
         "Справжня g1(x)", "Справжня g2(x)", 
         "Оцінка g1(x)", "Оцінка g2(x)"
       ), 
       col = c("green", "orange", "darkgreen", "yellow"), 
       lwd = line.thickness, 
       lty = c(2, 2, 1, 1))
par(mai=old.mai)

# Назва рисунку
error.types <- c("N(0, 0.0025)", "N(0, 1.25)", "T(5)")

title.name <- paste(
  ifelse(is.loc.lin.now, "LL", "NW"), ": ",
  "Приклад функції зв'язку №", example.num.func,
  ", розподіл похибок: ", error.types[example.num.err],
  sep=""
)

mtext(title.name, line=0, side=3, outer=TRUE, cex=1)

dev.off()
