from scipy.optimize import curve_fit
import numpy as np

#Fitting Model
#Assume that the input data is normalized and follows the model below:
def Model(x,gamma):
    global Beta
    return Beta*np.exp(-2*gamma*x)
#Import the normalized data to be fit by the model
def DLS_linear_fit(tau,g2,B,beta):
    global Beta
    Beta = beta
    g3 = g2
    #Lower baundry
    eps = 2.22*10**(-16)
    #Initial guess
    initialGuess = 1
    Params, Pcov = curve_fit(Model,tau,g3,p0 = initialGuess,bounds=((eps),(np.inf)))
    #Use the gamma to find the new data
    # Get data based on the fit
    Fit_data = Model(tau, Params[0])
    return Fit_data , Params
