
__author__= "Lucie Vigreux"

'''

This python code was developped during my internship at LARIS (http://laris.univ-angers.fr/fr/index.html). This code has been adapted  from the Matlab code ARFIMA Simulations of Simone Fatichi's (https://www.mathworks.com/matlabcentral/fileexchange/25611-arfima-simulations).

This code has only been tested to generate signals as described in Boris Podobnik and H. Eugene Stanley:  "Detrended Cross-Correlation Analysis: A New Method for Analyzing Two Nonstationary Time Series" (2008) (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.084102). Ie for signals with N fixed, 0<d<0.5 fixed, normal random noise : er =np.random.normal(0,1,N) and no other input.



This python code implements a function to generate ARFIMA (AutoRegressive Fractional Integral Moving Average) models. These models  generalize ARIMA (autoregressive integrated moving average) and ARMA (autoregressive moving average) models. ARFIMA models allow non-integer values of the differencing parameter and are useful in modeling time series with long memory.
The model is generally represented  as ARFIMA(p,d,q) model where d is the differencing parameter and p and q are the order of the autoregressive and moving average parts of the model respectively.

This package use the numpy package (http://numpy.scipy.org)

'''

import numpy as np
import math as m


# ARFIMA Simulation
#Inputs :

#N :length of the time series we want to generate
#F = [ F1 F2 F3 .... ] : Parameters of the AR model, length(F) is the order p. Default p = 0
#O = [ O1 O2 O3 .... ] : Parameters of the MA model, length(O) is the order q. Default q = 0
#d : fractionally differencing parameter, default d = 0
#stdx (optional input) : parameter to force the standard deviation of the output time series. Impose std(Z)==stdx
#er (optional input): predefined time seres of white noise

#Outputs :

#Z is the time series simulated with the ARFIMA model

def ARFIMA_SIM(N,F=[],O=[],d=0,stdx=np.nan,er=0):
    #inizialization
    X=np.zeros(N)
    Y=np.zeros(N)
    Z=np.zeros(N)
    F=np.array(F)
    O=np.array(O)
    e=np.random.normal(0,1,N)
    plt.plot(e)
    if type(er)!=int:
        e=er
    if (np.count_nonzero(O)==0) and (np.count_nonzero(F)==0) and (d==0) :
        Z=e
        return
    # N =length
    MA_ord=len(O)
    AR_ord=len(F)
    #Computing part: MA(q)
    t=0
    if MA_ord >= 1 :
        for t in range(N):
            j=0
            map=0
            for j in range(MA_ord):
                if t > j :
                    map+=O[j]*e[t-j]
            X[t]= e[t]+map
    else :
        X=e
    #Computing part: d
    if d==0 :
        Y=X
    else :
        infi=100
        s=0
        b=np.zeros(infi)
        for s in range(infi):
            b[s]=m.gamma(s-d)/(m.gamma(s+1)*m.gamma(-d))
        for t in range(N):
            Y[t]=0
            for s in range(infi):
                if t > s :
                    Y[t]+=b[s]*X[t-s]
    #Computing part: AR(p)
    if AR_ord >=1:
        for t in range(N):
            j=0
            arp=0
            for j in range(AR_ord):
                if t > j :
                    arp-=F[j]*Z[t-j]
            Z[t]=Y[t]+arp;
    else :
        Z=Y
    Z=np.transpose(Z)
    if not(np.isnan(stdx)) :
        Z=Z*stdx/np.std(Z)
    return Z
