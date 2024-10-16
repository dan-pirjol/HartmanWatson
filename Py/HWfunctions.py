import numpy as np
from numpy.polynomial.polynomial import polyval

# The series expansions for F(rho), G(rho) converge for rho in [0.03,32.88]
# Use tail expansions outside of this range
def Fexp(rho, n=5):
    # The coefficients in the expansions of F(rho) 
    log_rho = np.log(1/rho)
    cF = np.array([np.pi**2/2-1, 1., 1., -2/15, 19/525, -22/2625,
                   4742/3031875, -43636/197071875, 146287/6897515625])
    mid = polyval(log_rho, cF)
    
    low = 0.5*log_rho**2 + log_rho*np.log(np.abs(2*log_rho)) - log_rho + np.log(np.abs(2*log_rho))**2 + 0.5*np.pi**2
    high = rho + 0.5*np.pi**2/(1+rho)

    aux = np.where(rho < 0.03, low, mid)
    aux = np.where(rho > 32.0, high, aux)
    return aux

#----------------------------------------------------------

def Gexp(rho, n=5):
    # The coefficients in the expansions of G(rho) 
    log_rho = np.log(1/rho)
    cG = np.array([1., 1/5, -1/70, -1/1050, 299/323400, -96917/525525000, -107749/10032750000,
                   27333619/1876124250000, -308907281743/109790791110000000])
    mid = np.sqrt(3)*polyval(log_rho, cG)

    low = np.sqrt(np.abs(log_rho))*(1-rho**2) + 0.5*(np.log(np.abs(2*log_rho))+1)/np.sqrt(np.abs(log_rho))
    high = np.pi*rho/(1+rho)**(1.5)

    aux = np.where(rho < 0.03, low, mid)
    aux = np.where(rho > 32.0, high, aux)
    return aux

#-----------------------------------------------------------
# function returning theta(r,t)

def thetaHW(r, t, n=5):
    f1 = 1/(2*np.pi*t)
    rho = r*t
    aux = f1*Gexp(rho, n)*np.exp(-1/t*(Fexp(rho, n) - 0.5*np.pi**2))
    return aux
    