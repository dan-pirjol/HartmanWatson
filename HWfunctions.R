# The series expansions for F(rho), G(rho) converge for rho in [0.03,32.88]
# Use tail expansions outside of this range
Fexp <- function(rho,n=5){
  PI = 4*atan(1.0)
  # The coefficients in the expansions of F(rho) 
  cF <- c(1,1,-2/15,19/525,-22/2625,4742/3031875,-43636/197071875,146287/6897515625)
  
  mid <- PI^2/2-1
  for (i in 1:n) mid <- mid + cF[i]*log(1/rho)^i
  
  low <- 0.5*log(1/rho)^2 + log(1/rho)*log(abs(2*log(1/rho))) -log(1/rho) + log(abs(2*log(1/rho)))^2+0.5*PI^2
  high <- rho + 0.5*PI^2/(1+rho)
  
  if(rho < 0.03) aux = low
  else if (rho > 32.00) aux = high
  else aux = mid
  
  return(aux)
}

#----------------------------------------------------------
Gexp <- function(rho,n=5){
  PI = 4*atan(1.0)
  # The coefficients in the expansions of G(rho) 
  cG <- c(1/5,-1/70,-1/1050,299/323400,-96917/525525000,-107749/10032750000,27333619/1876124250000,
          -308907281743/109790791110000000)
  mid <- 1
  for (i in 1:n) mid <- mid + cG[i]*log(1/rho)^i
  mid <- sqrt(3)*mid
  
  low <- sqrt(abs(log(1/rho)))*(1-rho^2) + 0.5*(log(abs(2*log(1/rho)))+1)/sqrt(abs(log(1/rho)))
  high <- PI*rho/(1+rho)^(1.5)
  
  if(rho < 0.03) aux = low
  else if (rho > 32.00) aux = high
  else aux = mid
  
  return(aux)
}
#-----------------------------------------------------------
# function returning theta(r,t)

thetaHW <- function(r, t, n=5){
  f1 <- 1/(2*PI*t)
  rho <- r*t
  
  aux <- f1*Gexp(rho,n)*exp(-1/t*(Fexp(rho,n) - 0.5*PI^2))
  return(aux)
}
#-----------------------------------------------------------
# load pre-computed values of F,G on a grid [0.01 - 5.00] with step 0.01
allData <- read.csv(file = 'tableFGto5.csv')
#-----------------------------------------------------------
library("stats")
Finterp <- function(rho){
  z <- approx(allData$rho,allData$F,rho,method="linear",rule=2)

  f <- z$y
  
  return(f)
}
#-----------------------------------------------------------
Ginterp <- function(rho){
  z <- approx(allData$rho,allData$G,rho,method="linear",rule=2)
  
  g <- z$y
  
  return(g)
}
#-----------------------------------------------------------
