source("HWfunctions.R")

#--------------------------------------------------------
# plot the function F(rho)
x <- seq(0, 5, 0.01)
n <- length(x)

d <- c()
for (i in 1:n) d <- c(d,Ffunc(x[i]))

plot(d~x, type="l", col="blue", main="F(rho)", ylim=c(0,12))


#------------------------------------------------------
# Plot the function G(rho)
d <- c()
for (i in 1:n) d <- c(d,Gfunc(x[i]))

plot(d~x, type="l", col="blue", main="G(rho)")
#-------------------------------------------

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

head(y1)

