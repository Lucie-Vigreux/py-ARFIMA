# py-ARFIMA
This python code was developped during my internship at LARIS (http://laris.univ-angers.fr/fr/index.html). This code has been adapted  from the Matlab code ARFIMA Simulations of Simone Fatichi's (https://www.mathworks.com/matlabcentral/fileexchange/25611-arfima-simulations).

This code has only been tested to generate signals as described in Boris Podobnik and H. Eugene Stanley:  "Detrended Cross-Correlation Analysis: A New Method for Analyzing Two Nonstationary Time Series" (2008) (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.084102). Ie for signals with N fixed, 0<d<0.5 fixed, normal random noise : er =np.random.normal(0,1,N) and no other input.



This python code implements a function to generate ARFIMA (AutoRegressive Fractional Integral Moving Average) models. These models  generalize ARIMA (autoregressive integrated moving average) and ARMA (autoregressive moving average) models. ARFIMA models allow non-integer values of the differencing parameter and are useful in modeling time series with long memory.
The model is generally represented  as ARFIMA(p,d,q) model where d is the differencing parameter and p and q are the order of the autoregressive and moving average parts of the model respectively.

This package use the numpy package (http://numpy.scipy.org)
