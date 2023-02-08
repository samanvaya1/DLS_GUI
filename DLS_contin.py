import numpy as np
from scipy import optimize
import math

#Use the Contin Method

#Function to truncate a number to get a faster algorithuim run time
#Functions for Contin
def integral_con(x):
    c = []
    ceq = np.sum(x)-1
    return ceq

def msd(a,g,A,alpha,om,gamma,W):
    E = ((np.sum(np.square(np.subtract((np.dot(A,a)),g)))))
    R = np.sum(np.square(np.subtract(np.add(a[0:len(a)-3],a[2:len(a)-1]),np.multiply(2,a[1:len(a)-2]))))
    return E+(alpha**2)*R

def DLS_contin(tau,g1,gamma,alpha,W,x0):
    #Vector for initial guess for forier coefficients
    #Vector to be used as initial guess for Fourier coefficients
    #Findning A, the transformation matrix defined by exp(-tM*gammaN)
    [gammaN , tM] = np.meshgrid(gamma,tau)
    A = np.exp(-2*np.multiply(tM,gammaN))
    #U,s,V = np.linalg.svd(A)
    Row,Col = gammaN.shape
    #s_holder = np.zeros((Row,Col))
    #np.fill_diagonal(s_holder,s)
    #s = s_holder
    #Cut off for the eigenvalue and apply to matrix
    #eigenvalue_cutoff = 0
    #k = 0
    #while(k<Row and k <Col):
     #   if(s[k,k]<eigenvalue_cutoff):
      #      s[k,k] =0
       # k = k+1
    #Compute the optimal solution to truncated SVD problem
    #A = np.dot(np.dot(U,s),V)
    omCol = x0.shape
    om = np.zeros(omCol)

    #Constraints for minimization
    x0 = x0.transpose()
    lb_fmin = np.zeros(len(x0))
    ub_fmin = np.inf*np.ones(len(x0))
    non_lin  = lambda x: integral_con(x)

    #The nonlinear constraint for the function
    con = {'type': 'eq', 'fun': non_lin}

    A_fmin = []
    b_fmin = []
    alpha = float(alpha)
    test_fn = lambda x: msd(x,g1,A,alpha,om,gamma,W)
    bnds = optimize.Bounds(lb_fmin,ub_fmin,keep_feasible=True)
    x_opt = optimize.minimize(test_fn,x0,method='SLSQP',bounds = bnds,constraints=con ,options={'ftol':10**-14,'maxiter':1000000})
    #TNC a bit slow bad size dist and meh error
    #SLSQP pretty good speed and stats
    x_min = x_opt.x
    fval = test_fn(x_min)
    fit_data = np.dot(A,x_min)
    return fit_data , gamma , x_min