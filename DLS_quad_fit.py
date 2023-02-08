import numpy as np
from scipy.optimize import curve_fit

#Assume that the input data is normalized and follows the model below:
def Model(x,gamma,k2):
    global Beta
    return Beta*np.exp(-2*gamma*x)*np.square(1+0.5*np.square(k2*x))

#The function to do the numerical analysis
def DLS_quad_fit(tau,g2,B,beta,initialGuess):
    #Globalize beta because it is a constant
    global Beta
    g3 = g2
    #Get Beta as the first two points
    Beta = beta
    #Lower baundry
    eps = 2.22*10**(-16)

    Params, Pcov = curve_fit(Model,tau,g3,p0 = initialGuess,bounds=((eps,eps),(np.inf,np.inf)))
    #Get data based on the fit
    Fit_data = Model(tau,Params[0],Params[1])
    return Params, Fit_data