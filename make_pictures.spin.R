library(spatstat)

data(cells)
plot(cells, main = "")

set.seed(1234)
csr <- rpoispp(cells$n, win = cells$window)
plot(csr)

op <- par(mfrow = c(1,1))
par(mfrow = c(1,2),
    mar = c(0.1,0.1,0.1,0.1))
plot(cells, main = "")
plot(rpoispp(cells$n, win = cells$window), main = "")
par(op)


# Create covariate gradient
grad_lambda <- function(x,y) 100*x 
grad_pp <- rpoispp(grad_lambda, win = cells$window)
plot(grad_pp)

xs <- seq(0, 1, length.out = 200)
ys <- seq(0, 1, length.out = 200)
df <- as.data.frame(expand.grid(x = xs, y = ys))
df$lambda <- grad_lambda(xs,ys)

library(ggplot2)
ggplot(df) + 
  geom_tile(aes(x=x,y=y, fill = lambda))  

grad_lambda <- function(x,y) 100*x 
grad_pp <- rpoispp(grad_lambda, win = cells$window)
pp_df <- data.frame(x = grad_pp$x, y = grad_pp$y)

xs <- seq(0, 1, length.out = 200)
ys <- seq(0, 1, length.out = 200)
df <- as.data.frame(expand.grid(x = xs, y = ys))
df$lambda <- grad_lambda(xs,ys)
ggplot() + 
  geom_tile(data = df, aes(x=x,y=y, fill = lambda)) +
  geom_point(data = pp_df, aes(x=x, y=y)) +
  theme_void() +
  guides(fill = FALSE)

hill_lambda <- function(x,y) (x - 0.5)^2 + (y - 0.5)^2
lambda2 <- function(x,y) grad_lambda(x,y) + 300*hill_lambda(x,y)

pp2 <- rpoispp(lambda2, win = cells$window)
pp_df <- data.frame(x = pp2$x, y = pp2$y)

xs <- seq(0, 1, length.out = 200)
ys <- seq(0, 1, length.out = 200)
df <- as.data.frame(expand.grid(x = xs, y = ys))
df$lambda <- lambda2(xs,ys)
ggplot() + 
  geom_tile(data = df, aes(x=x,y=y, fill = lambda)) +
  geom_point(data = pp_df, aes(x=x, y=y)) +
  theme_void() +
  guides(fill = FALSE)

set.seed(100)
lgcp_pp <- rLGCP(model = "exp", 
                 mu = function(x,y) log(grad_lambda(x,y)), 
                 win = cells$window,
                 param = list(var = 0.2, scale = 0.1))
lgcp_pp_df <- data.frame(x = lgcp_pp$x, y = lgcp_pp$y)
xs <- seq(0, 1, length.out = 200)
ys <- seq(0, 1, length.out = 200)
df <- as.data.frame(expand.grid(x = xs, y = ys))
df$lambda <- grad_lambda(xs,ys)
ggplot() + 
  geom_tile(data = df, aes(x=x,y=y, fill = lambda)) +
  geom_point(data = lgcp_pp_df, aes(x=x, y=y)) +
  theme_void() +
  guides(fill = FALSE)


library(RandomFields)
## first let us look at the list of implemented models
RFgetModelNames(type="positive definite", domain="single variable",
                iso="isotropic")
## our choice is the exponential model;
## the model includes nugget effect and the mean:
model <- RMexp(var=2, scale=10) + RMtrend(mean=0)
## define the locations:
from <- 0
to <- 20
x.seq <- seq(from, to, length=200)
simu <- RFsimulate(model, x=x.seq)
plot(simu, ylim = c(-3,3))


class(simu)
str(simu)

library(ggplot2)
df <- data.frame(x = x.seq, y = simu@data$variable1)
ggplot(df) +
  geom_line(aes(x=x, y =y)) +
  ylim(c(-5,5)) +
  xlab("s") + 
  ylab("f")


df$y
df$x




# Exp correlation fn
cexp <- function(r, scale) exp(-1/(2*scale)*r)

sc10 <- 10
df10 <- data.frame(x = x.seq, cexp = cexp(x.seq, sc10))

sc4 <- 4
df4 <- data.frame(x = x.seq, cexp = cexp(x.seq, sc4))

sc1 <- 1
df1 <- data.frame(x = x.seq, cexp = cexp(x.seq, sc1))

ggplot() +
  geom_line(data = df10, aes(x=x, y=cexp, colour = "black")) +
  geom_line(data = df4, aes(x=x, y=cexp, colour = "orange")) +
  geom_line(data = df1, aes(x=x, y=cexp, colour = "blue")) +
  scale_colour_manual(values = c("black", "red", "blue"),
                      labels = c("scale = 10", "scale = 4", "scale = 1")) +
  xlab("r") +
  ylab(expression(paste(kappa, "(r)")))


ggplot() +
  geom_line(data = df10, aes(x=x, y=cexp)) +
  geom_line(data = df1, aes(x=x, y=cexp))

ggplot() +
  geom_line(data = df10, aes(x=x, y=cexp))



# Simulate now from these
m10 <- RMexp(var=2, scale=sc10) + RMtrend(mean=0)
m4 <- RMexp(var=2, scale=sc4) + RMtrend(mean=0)
m1 <- RMexp(var=2, scale=sc1) + RMtrend(mean=0)

sim10 <- RFsimulate(m10, x=x.seq)
sim4 <- RFsimulate(m4, x=x.seq)
sim1 <- RFsimulate(m1, x=x.seq)

dfs10 <- data.frame(x = x.seq, y = sim10@data$variable1)
dfs4 <- data.frame(x = x.seq, y = sim4@data$variable1)
dfs1 <- data.frame(x = x.seq, y = sim1@data$variable1)
ggplot() +
  geom_line(data = dfs10, aes(x=x, y=y, colour = "black")) +
  geom_line(data = dfs4, aes(x=x, y=y, colour = "orange")) +
  geom_line(data = dfs1, aes(x=x, y=y, colour = "blue")) +
  scale_colour_manual(values = c("black", "red", "blue"),
                      labels = c("scale = 10", "scale = 1", "scale = 4")) +
  xlab("s") +
  ylab("f")
