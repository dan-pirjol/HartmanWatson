source("HWfunctions.R")

#--------------------------------------------------------
# plot the function F(rho)
x <- seq(0, 5, 0.01)
n <- length(x)

d <- c()
for (i in 1:n) d <- c(d,Fexp(x[i]))

plot(d~x, type="l", col="blue", main="F(rho)", ylim=c(0,12))


#------------------------------------------------------
# Plot the function G(rho)
d <- c()
for (i in 1:n) d <- c(d,Gexp(x[i]))

plot(d~x, type="l", col="blue", main="G(rho)")
#------------------------------------------------------

# plot the function theta(r,t)
# Reproduce Figure 4 in the paper
x <- seq(0, 4, 0.01)
n <- length(x)

y1 <- c()
y2 <- c()
y3 <- c()
for (i in 1:n) {
  y1 <- c(y1,thetaHW(0.5,x[i]))
  y2 <- c(y2,thetaHW(1.0,x[i]))
  y3 <- c(y3,thetaHW(1.5,x[i]))
}

plot(y1~x, type="l", col="black", main="thetaHW(r,t)", ylim=c(0,2.5))
lines(y2~x, type="l", col="blue")
lines(y3~x, type="l", col="red")

#---------------------------------------------------------
# compare the series expansion with the exact evaluations

x <- seq(0.5, 2.0, 0.01)
n <- length(x)

f1 <- c()
f2 <- c()

for (k in 1:n) {
  xarg <- x[k]
  f1 <- c(f1, Finterp(xarg))
  f2 <- c(f2, Fexp(xarg,5))
}
                     
f3 <- f1 - f2

#fd <- data.frame("rho"=x, "exact"=f1, "series"=f2, "diff"=f3, "log"=log(abs(f3),10))
#fd

plot(log(abs(f3),10) ~ x, type="l", main="F(rho) log-error (n=5,6,7)", xlab="rho",ylim=c(-16,0))
lines(log(abs(f3),10) ~ x, type="l", col="blue")
lines(log(abs(f3),10) ~ x, type="l", col="red")

#----------------------------------------------------------
g1 <- c()
g2 <- c()

for (k in 1:n) {
  xarg <- x[k]
  g1 <- c(g1, Ginterp(xarg))
  g2 <- c(g2, Gexp(xarg,5))
}

g3 <- g2 - g1

gd <- data.frame("rho"=x, "exact"=g1, "series"=g2, "diff"=g3)
gd

plot(log(abs(g3),10) ~ x, type="l", main="G(rho) log-error", xlab="rho")
lines(log(abs(g3),10) ~ x, type="l", col="blue")
lines(log(abs(g3),10) ~ x, type="l", col="red")

#----------------------------------------------------------
