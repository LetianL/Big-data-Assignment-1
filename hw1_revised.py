# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
from math import *
from scipy.stats import norm

def d1(S,K,T,r,sigma):
    return(log(S/K)+(r+sigma**2/2.)*T)/(sigma*sqrt(T))
def d2(S,K,T,r,sigma):
    return d1(S,K,T,r,sigma)-sigma*sqrt(T)

def bs_call(S,K,T,r,sigma):
    return S*norm.cdf(d1(S,K,T,r,sigma))-K*exp(-r*T)*norm.cdf(d2(S,K,T,r,sigma))
  
def bs_put(S,K,T,r,sigma):
    return K*exp(-r*T)-S+bs_call(S,K,T,r,sigma)

def delta_call(S,K,T,r,sigma):
    return norm.cdf(d1(S,K,T,r,sigma))
def Gamma_call(S,K,T,r,sigma):
    return norm.pdf(d1(S,K,T,r,sigma))/(S*sigma*np.sqrt(T))
def Vega_call(S,K,T,r,sigma):
    return S*norm.pdf(d1(S,K,T,r,sigma))*np.sqrt(T)
def Theta_call(S,K,T,r,sigma):
    aux1 = -S*norm.pdf(d1(S,K,T,r,sigma))*sigma/(2*np.sqrt(T))
    aux2 = -r*K*np.exp(-r*T)*norm.cdf(d2(S,K,T,r,sigma))
    return aux1+aux2
def Rho_call(S,K,T,r,sigma):
    return K*T*np.exp(-r*T)*norm.cdf(d2(S,K,T,r,sigma))

def delta_put(S,K,T,r,sigma):
    return delta_call(S,K,T,r,sigma) - 1
def Gamma_put(S,K,T,r,sigma):
    return norm.pdf(d1(S,K,T,r,sigma))/(S*sigma*np.sqrt(T))
def Vega_put(S,K,T,r,sigma):
    return S*norm.pdf(d1(S,K,T,r,sigma)*np.sqrt(T))
def Theta_put(S,K,T,r,sigma):
    aux1 = -S*norm.pdf(d1(S,K,T,r,sigma))*sigma/(2*np.sqrt(T))
    aux2 = r*K*np.exp(-r*T)*norm.cdf(-d2(S,K,T,r,sigma))
    return aux1+aux2
def Rho_put(S,K,T,r,sigma):
    return -K*T*np.exp(-r*T)*norm.cdf(-d2(S,K,T,r,sigma))

if __name__ == "__main__":
    print('BS_Call is: '+ str(bs_call(100,100,1,0,0.3))+'\n')
    print('BS_put is: '+str(bs_put(100,100,1,0,0.3))+'\n')  
    print('delta_call is: '+str(delta_call(100,100,1,0,0.3))+'\n') 
    print('Gamma_call is: ' + str(Gamma_call(100,100,1,0,0.3))+'\n')
    print('Vega_call is: '+ str(Vega_call(100,100,1,0,0.3))+'\n')
    print('Theta_call is: '+ str(Theta_call(100,100,1,0,0.3))+'\n')
    print('Rho_call is: '+str(Rho_call(100,100,1,0,0.3))+'\n')
    print('\n')
    print('delta_put is: '+str(delta_put(100,100,1,0,0.3))+'\n') 
    print('Gamma_put is: ' + str(Gamma_put(100,100,1,0,0.3))+'\n')
    print('Vega_put is: '+ str(Vega_put(100,100,1,0,0.3))+'\n')
    print('Theta_put is: '+ str(Theta_put(100,100,1,0,0.3))+'\n')
    print('Rho_put is: '+ str(Rho_put(100,100,1,0,0.3))+'\n')



  
      
     
