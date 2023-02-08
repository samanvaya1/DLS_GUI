import DLS_cubic_fit
import DLS_quad_fit
import DLS_linear_fit
import numpy as np
from scipy.optimize import curve_fit

#Assume that the input data is normalized and follows the model below:
def Model(x,gamma,beta,k2,k3,k4):
    return beta*np.exp(-2*gamma*x)*np.square(1+0.5*np.square(k2*x)+(1/6)*np.power(k3*x,3)+(1/24)*np.power(k4*x,4))

#The function to do the numerical analysis
def DLS_quartic_fit(tau,g2):
    #Globalize beta because it is a constant
    beta_guess = 0.5*(g2[0]+g2[1])
    Beta = beta_guess
    g3 = g2
    B =0.1
    #Lower baundry
    eps = 2.22*10**(-16)
    #Initial guess
    # Initial guess
    Params1, Fit = DLS_linear_fit.DLS_linear_fit(tau, g3, B, Beta)
    initialGuess = [Params1[0], 1]
    OldParams, OldData = DLS_quad_fit.DLS_quad_fit(tau, g3, B, Beta, initialGuess)
    initialGuess = [Params1[0],Params1[1], 1]
    OldParams, OldData = DLS_cubic_fit.DLS_cubic_fit(tau, g3, B, Beta,initialGuess)
    initialGuess = [OldParams[0],beta_guess,OldParams[1],OldParams[2],1]
    Params, Pcov = curve_fit(Model,tau,g3,p0 = initialGuess,bounds=((eps,eps,eps,eps,eps),(np.inf,np.inf,np.inf,np.inf,np.inf)))

    #Get data based on the fit
    Fit_data = Model(tau,Params[0],Params[1],Params[2],Params[3],Params[4])
    return Params, Fit_data