import numpy as np
import scipy.stats

def f(theta, x, returngrad=False):
    fval = (theta[0])*(theta[1] - x)
    if returngrad==False:
        return fval
    else:
        if (np.isscalar(x)):
            xshape = ()
        else:
            xshape = x.shape
        dfdtheta = np.zeros(theta.shape + xshape)
        # want this to broadcast correctly regardless of whether x is a matrix or vector
        dfdtheta[0,...] = fval/theta[0]
        dfdtheta[1,...] = theta[0]
        dfdtheta[2,...] = 0
        return fval, dfdtheta
    

def g(theta, x, returngrad=False):
    if (np.isscalar(x)):
        xshape = ()
        gval = theta[2]
    else:
        xshape = x.shape
        gval = np.ones(xshape)*theta[2]
    if returngrad==False:
        return gval
    else:
        dgdtheta = np.zeros(theta.shape + xshape)
        # want this to broadcast correctly regardless of whether x is a matrix or vector
        dgdtheta[0,...] = 0
        dgdtheta[1,...] = 0
        dgdtheta[2,...] = 1


def integrandmat(x, y, h, theta):
    fval = f(theta,y)
    gval = np.abs(g(theta,y))
    mu = y + fval*h
    sd = gval*np.sqrt(h)
    out = scipy.stats.norm.pdf(x, loc=mu, scale=sd)
    return(out)


def dtqall(theta, h, k, M, numsteps, init, final):
    grid = np.arange(-M,(M+1))*k
    gridx, gridy = np.meshgrid(grid, grid, sparse=False, indexing='ij')
    A = integrandmat(gridx, gridy, h, theta)
    
    # supposed to be a column vector
    lamb = np.array([integrandmat(grid, init, h, theta)]).T
    
    # supposed to be a row vector
    gamm = np.array([k * integrandmat(final, grid, h, theta)])
    
    lambdalist = []
    lambdalist.append(lamb)
    gammalist = []
    gammalist.append(gamm)
    
    for j in range(numsteps-2):
        lambdalist.append(np.dot(k*A,lambdalist[j]))
        gammalist.append(np.dot(gammalist[j],k*A))
    
    complete = np.asscalar(np.dot(gammalist[0],lambdalist[numsteps-2]))
    
    # both first and list wil be column vectors
    first = (1/k)*((gammalist[numsteps-2].T) * lambdalist[0])/complete
    last = (1/k)*((gammalist[0].T) * lambdalist[numsteps-2])/complete
    
    # compute all intermediate 2d pdfs
    pdf2dlist = []
    for j in range(numsteps-2):
        thisone = ((gammalist[(numsteps-3-j)] * lambdalist[j]).T * A)/(k*complete)
        pdf2dlist.append(thisone)

    return complete, first, last, pdf2dlist


# test run
# k = 0.01
# M = 800
# deltat = 1.
# numsteps = 50
# h = deltat/numsteps
# theta = np.array([1., 0., 2.])
# init = 1.
# final = 3.
# complete = dtqall(theta, h, k, M, numsteps, init, final)


# write all matrices to disk to compare with R code
# import scipy.io
# scipy.io.mmwrite(a=complete[1],target="first.mtx")
# scipy.io.mmwrite(a=complete[2],target="last.mtx")
# for j in range(len(complete[3])):
#     fname = "pdf2d"+str(j)+".mtx"
#     scipy.io.mmwrite(a=complete[3][j],target=fname)




