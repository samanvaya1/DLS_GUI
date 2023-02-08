import numpy as np
from scipy import optimize


def msd(x,g1,A,alpha):
    x = np.abs(x)
    VAR = (np.sum(np.square(np.subtract(g1,np.dot(A,x)))))
    #Om_ega = np.diff(np.diff(x))
    #om = np.zeros(len(x)-2)
    #REG = (alpha**2) * (np.sum(np.square(np.subtract(om,Om_ega))))
    REG = 0
    #(np.square(x[0]) + np.square(2*x[0]-x[1]) + np.sum(np.square(np.subtract(np.add(x[0:len(a)-4],x[2:len(x)-2]),np.multiply(2,x[1:len(x)-3])))) + np.square(np.multiply(x[len(x)-1],-2)+ x[len(x)-2])+np.square(x[len(x)-1]))

    return VAR + REG

def integral_con(x):
    c = []
    ceq = np.sum(x)-1
    return ceq

def DLS_DYNALS(t,g,gamma,x0):

    g1 = g
    alpha =0.5

    #Find transformation matrix
    [gammaN , tM] = np.meshgrid(gamma,t)
    A = np.exp(-np.multiply(tM,gammaN))
    U,s,V = np.linalg.svd(A)
    Row,Col = gammaN.shape
    s_holder = np.zeros((Row,Col))
    np.fill_diagonal(s_holder,s)
    s = s_holder

    #Cut off for the eigenvalue and apply to matrix
    eigenvalue_cutoff = 0.2
    k = 0
    while(k<Row and k <Col):
        if(s[k,k]<eigenvalue_cutoff):
            s[k,k] = 0
        k = k +1
    #FInd s pseudo
    s_pseudo_inv =s
    k =0
    while(k<Row and k <Col):
        if(s[k,k]==0):
            s_pseudo_inv =0
        else:
            s_pseudo_inv = 1/(s[k,k])
        k = k+1

    A = np.dot(np.dot(U,s),V)


    #Constraints for minimization
    x0 = x0.transpose()
    lb_fmin = np.zeros(len(x0))
    ub_fmin = np.inf*np.ones(len(x0))
    non_lin  = lambda x: integral_con(x)

    #The nonlinear constraint for the function
    con = {'type': 'eq', 'fun': non_lin}

    A_fmin = []
    b_fmin = []

    test_fn = lambda x: msd(x,g1,A,alpha)
    bnds = optimize.Bounds(lb_fmin, ub_fmin, keep_feasible=True)
    x_opt = optimize.minimize(test_fn, x0, method='SLSQP', bounds=bnds, constraints=con,options={'ftol': 10 ** -14, 'maxiter': 10000, 'disp': False})
    x_min = x_opt.x
    fval = test_fn(x_min)
    fit_data = np.dot(A, x_min)
    x = x_min / np.sum(x_min)


    #Adjust for the artifact tail end of the fit. Typicall occurs within the last two points but will do last 5 percent of points just to be sure the residual
    NumPoints = int(0.05*len(x))
    for index in range(NumPoints):
        x[len(x) - 1 - index] = 10**-7

    x = x / np.sum(x)
    return fit_data,gamma,x
