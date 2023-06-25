# HartmanWatson
R functions for evaluation of the Hartman-Watson function $$\theta(r,t) = \frac{r}{\sqrt{2\pi^3 t}} e^{\frac{\pi^2}{2t}}\int_0^\infty e^{-\frac{\xi^2}{2t}} e^{-r\cosh \xi} \sinh \xi \sin \frac{\pi \xi}{t} d\xi$$

The function $\theta(r,t)$ is evaluated as the leading order term in the $t\to 0$ asymptotic expansion at fixed $r t = \rho$ proposed in [Pirjol (2020)](https://arxiv.org/abs/2001.09579).
The leading term in this expansion is 
$$\theta(\rho/t,t)= \frac{1}{2\pi t} e^{-\frac{1}{t}(F(\rho) - \frac{\pi^2}{2}} G(\rho)(1+O(t))$$ 

The functions $F(\rho),G(\rho)$ appearing in this expansion are known exactly. The $O(t)$ term is also known in closed form.

Denoting the leading term in this expansion $\hat \theta(\rho/t,t)$, the $O(t)$ error is bounded as
$|\theta(\rho/t,t) - \hat\theta(\rho/t,t)| \leq \frac{1}{70} t \hat \theta(\rho/t,t)$ uniformly over $\rho$.

The code approximates the functions $F(\rho),G(\rho)$ as series in $\log(1/\rho)$ using the methods described in [Nandori, Pirjol (2021)](https://arxiv.org/abs/2209.09412).
These series converge within the convergence domain $|\log\rho| < 3.49295$. Outside of this region, the tail asymptotics of $F(\rho),G(\rho)$ are used. 

The function **Fexp(rho,n)** computes the function $F(\rho)$ using an expansion in $\log(1/\rho)$ keeping $n$ terms (default $n=5$, maximum allowed value 8). 

The function **Gexp(rho,n)** computes the function $G(\rho)$ using an expansion in $\log(1/\rho)$ keeping $n$ terms (default $n=5$, maximum allowed value 8).

The function **thetaHW(r,t,n)** returns $\hat \theta(r,t)$ using the expansions for $F(\rho),G(\rho)$.

The function **Finterp(rho)** evaluates the function $F(\rho)$ by linear interpolation from a table of pre-computed exact values on a grid $\rho:[0.01,5.00]$ with step 0.01. Flat extrapolation outside the grid range.

The function **Ginterp(rho)** evaluates $G(\rho)$ by linear interpolation from a table of pre-computed exact values on a grid $\rho:[0.01,5.00]$ with step 0.01. Flat extrapolation outside the grid range.
 

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



