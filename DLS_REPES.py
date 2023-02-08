import numpy as np
from scipy import optimize


def error(A,x,g1,alpha):
    x = np.abs(x)
    VAR = np.sum(np.square(np.subtract(g1,np.dot(A,x))))

    #Om_ega = np.subtract(np.add(x[0:len(x)-3],x[2:len(x)-1]),np.multiply(2,x[1:len(x)-2]))
    #Om_ega = (np.diff(x))
    #om = np.zeros(len(x)-1)
    #REG = (np.sum(np.square(om-Om_ega)))

    REG = (np.square(x[0]) + np.square(2*x[0]-x[1]) + np.sum(np.square(np.subtract(np.add(x[0:len(x)-4],x[2:len(x)-2]),np.multiply(2,x[1:len(x)-3])))) + np.square(np.multiply(x[len(x)-1],-2)+ x[len(x)-2])+np.square(x[len(x)-1]))

    return VAR+(alpha**2)*REG

def DLS_REPES(t,g,Tau,reg_param,a0):
    alpha = reg_param
    [Tau_M,tM] = np.meshgrid(Tau,t)
    A = np.exp(-tM/Tau_M)

    #Constraints for minimization
    lb_fmin = np.zeros(len(Tau))
    ub_fmin = np.inf*np.ones(len(Tau))
    A_fmin = []
    b_fmin = []

    test_fn = lambda x: error(A,x,g,alpha)
    bnds = optimize.Bounds(lb_fmin, ub_fmin, keep_feasible=True)
    a_opt = optimize.minimize(test_fn, a0, method='SLSQP', bounds=bnds,options={'ftol': 10 ** -14, 'maxiter': 1000000, 'disp': False})
    a = a_opt.x
    g_calc = np.dot(A,a)
    fit_data = g_calc
    a = np.divide(a,np.sum(a))
    # Adjust for the artifact tail end of the fit. Typicall occurs within the last two points but will do last 5 percent of points just to be sure the residual
    NumPoints = int(0.05 * len(a))

    for index in range(NumPoints):
        # Set Condition for when it is to be considered an artifact
        a[len(a) - 1 - index] = 10 ** -7

    a = a / np.sum(a)
    return fit_data,Tau,a
