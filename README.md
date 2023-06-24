# HartmanWatson
R functions for evaluation of the Hartman-Watson function $$\theta(r,t) = \frac{r}{\sqrt{2\pi^3 t}} e^{\frac{\pi^2}{2t}}\int_0^\infty e^{-\frac{\xi^2}{2t}} e^{-r\cosh \xi} \sinh \xi \sin \frac{\pi \xi}{t} d\xi$$

The function $\theta(r,t)$ is evaluated as the leading order term in the asymptotic expansion in [Pirjol (2020)](https://arxiv.org/abs/2001.09579).
$$\theta(\rho/t,t)= \frac{1}{2\pi t} e^{-\frac{1}{t}(F(\rho) - \frac{\pi^2}{2}} G(\rho)(1+O(t))$$ 

The functions $F(\rho),G(\rho)$ appearing in this expansion are approximated as series in $\log(1/\rho)$ using the methods described in [Nandori, Pirjol (2021)](https://arxiv.org/abs/2209.09412).
These series converge within the convergence domain $|\log\rho| < 3.49295$. Outside of this region, the tail asymptotics of $F(\rho),G(\rho)$ are used. 

The function **Ffunc(rho,n)** computes the function $F(\rho)$ using an expansion in $\log(1/\rho)$ keeping $n$ terms (default $n=5$, maximum allowed value 8). 

The function **Gfunc(rho,n)** computes the function $G(\rho)$ using an expansion in $\log(1/\rho)$ keeping $n$ terms (default $n=5$, maximum allowed value 8).

The function **thetaHW(r,t)** computes $\theta(r,t)$ using the expansions for $F(\rho),G(\rho)$ keeping $n=5$ terms.
This function returns $\hat \theta(r,t)$, the leading term in the $t\to 0$ asymptotic expansion of $\theta(r,t)$ at fixed $r t = \rho$. The error of this approximation is bounded as
$|\theta(r,t) - \hat\theta(r,t)| \leq \frac{1}{70} t \hat \theta(r,t)$ 

## **Sample usage**
```
#Reproduce the plots in Figure 4 of https://arxiv.org/abs/2001.09579
source("HWfunctions.R")
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

plot(y1~x, type="l", col="black", main="thetaHW(r,t)", ylim=c(0,2.5), xlab="t",ylab="theta(r,t)")
lines(y2~x, type="l", col="blue")
lines(y3~x, type="l", col="red")
```
![figure4](https://github.com/dan-pirjol/HartmanWatson/assets/60016102/7158ad3e-0d8d-411f-82ea-e8f3673b6a34)



