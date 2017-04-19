import sys
import ctypes
# import logging
import multiprocessing as mp

import numpy as np
import scipy.stats
import scipy.optimize

from contextlib import closing

# info = mp.get_logger().info

class gl:
    k = 0.05
    M = 100
    deltat = 0.1
    numsteps = 10
    h = deltat/numsteps
    grid = np.arange(-M,(M+1))*k
    size = 2*M + 1
    size2 = size*size
    gridx = mp.Array(ctypes.c_double, size2)
    gridy = mp.Array(ctypes.c_double, size2)


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
        return gval, dgdtheta


# using shared grids
def integrandmat(theta, arr):
    gridxvec = tonumpyarray(gl.gridx)
    gridx = gridxvec.view()
    gridx.shape = (gl.size, gl.size)
    gridyvec = tonumpyarray(gl.gridy)
    gridy = gridyvec.view()
    gridy.shape = (gl.size, gl.size)

    fval = f(theta,gridy)
    gval = np.abs(g(theta,gridy))
    mu = gridy + fval*gl.h
    sd = gval*np.sqrt(gl.h)
    arr[:] = scipy.stats.norm.pdf(gridx, loc=mu, scale=sd)


def integrandmat0(x, y, h, theta):
    fval = f(theta,y)
    gval = np.abs(g(theta,y))
    mu = y + fval*h
    sd = gval*np.sqrt(h)
    return scipy.stats.norm.pdf(x, loc=mu, scale=sd)


def dtqall(x):
    init = x[0]
    final = x[1]
    theta = x[2:5]

    # grid = np.arange(-M,(M+1))*k
    # gridx, gridy = np.meshgrid(grid, grid, sparse=False, indexing='ij')
    # A = integrandmat(gridx, gridy, h, theta)
    Avec = tonumpyarray(shared_prop)
    A = Avec.view()
    A.shape = (gl.size, gl.size)
    
    # supposed to be a column vector
    lamb = np.array([integrandmat0(gl.grid, init, gl.h, theta)]).T
    
    # supposed to be a row vector
    gamm = np.array([gl.k * integrandmat0(final, gl.grid, gl.h, theta)])
    
    lambdalist = []
    lambdalist.append(lamb)
    gammalist = []
    gammalist.append(gamm)
    
    for j in range(gl.numsteps-2):
        lambdalist.append(np.dot(gl.k*A,lambdalist[j]))
        gammalist.append(np.dot(gammalist[j],gl.k*A))
    
    complete = np.asscalar(np.dot(gammalist[0],lambdalist[gl.numsteps-2]))
    
    # both first and list wil be column vectors
    first = (1/gl.k)*((gammalist[gl.numsteps-2].T) * lambdalist[0])/complete
    last = (1/gl.k)*((gammalist[0].T) * lambdalist[gl.numsteps-2])/complete
    
    # compute all intermediate 2d pdfs
    pdf2dlist = []
    for j in range(gl.numsteps-2):
        thisone = ((gammalist[(gl.numsteps-3-j)] * lambdalist[j]).T * A)/(gl.k*complete)
        pdf2dlist.append(thisone)

    return complete, first, last, pdf2dlist


# grid version only
def logG(theta, lg, gradlg):
    gridxvec = tonumpyarray(gl.gridx)
    x = gridxvec.view()
    x.shape = (gl.size, gl.size)
    gridyvec = tonumpyarray(gl.gridy)
    y = gridyvec.view()
    y.shape = (gl.size, gl.size)

    fv, fgrad = f(theta,y,returngrad=True)
    gv, ggrad = g(theta,y,returngrad=True)
    pr = x - y - fv*gl.h
    pr2 = np.square(pr)
    gv2 = np.square(gv)
    lg[:] = -0.5*np.log(2*np.pi*gv2*gl.h) - pr2/(2*gv2*gl.h)

    # gradient of logG w.r.t. theta
    dlgdf = pr/gv2
    dlgdg = -1./gv + pr2/(gv2*gv*gl.h)
    
    n = len(theta)
    for i in range(n):
        gradlg[i,...] = dlgdf*fgrad[i,...] + dlgdg*ggrad[i,...]
        

def logG0(x, y, theta, h):
    fv, fgrad = f(theta,y,returngrad=True)
    gv, ggrad = g(theta,y,returngrad=True)
    pr = x - y - fv*h
    pr2 = np.square(pr)
    gv2 = np.square(gv)
    lg = -0.5*np.log(2*np.pi*gv2*h) - pr2/(2*gv2*h)
    dlgdf = pr/gv2
    dlgdg = -1./gv + pr2/(gv2*gv*h)
    
    n = len(theta)
    if np.isscalar(x):
        # in this case, y had better be a 1-d array
        gradlg = np.zeros(theta.shape + y.shape)
    else:
        # in this case, x had better be a 1-d or 2-d array
        gradlg = np.zeros(theta.shape + x.shape)
        
    for i in range(n):
        gradlg[i,...] = dlgdf*fgrad[i,...] + dlgdg*ggrad[i,...]
        
    return lg, gradlg
    

# allout is the result of a call to dtqall
# here x is the data
def qfun(allout):
    init = allout[0][0]
    final = allout[0][1]
    theta = allout[0][2:5]
    q = 0
    n = len(theta)
    gradq = np.zeros(n)

    # first term in the summation (j=1 case)
    part1 = logG0(gl.grid,init,theta,gl.h)
    q += np.asscalar(np.dot(part1[0],allout[1][1]))*gl.k

    # loop over components of theta
    for i in range(n):
        gradq[i] += np.asscalar(np.dot(part1[1][i,:],allout[1][1]))*gl.k


    gradqtest = (np.dot(part1[1],allout[1][1])*gl.k).ravel()
    print(np.sum(np.abs(gradq - gradqtest)))

    # last term in the summation (j=F case)
    part2 = logG0(final,gl.grid,theta,gl.h)
    q += np.asscalar(np.dot(part2[0],allout[1][2]))*gl.k

    # loop over components of theta
    for i in range(n):
        gradq[i] += np.asscalar(np.dot(part2[1][i,:],allout[1][2]))*gl.k

    # all intermediate terms
    # part3 = logG(gridx,gridy,theta,gl.h)
    arr = tonumpyarray(shared_part3)
    part3 = arr.view()
    part3.shape = (gl.size, gl.size)
    arrgrad = tonumpyarray(shared_part3grad)
    part3grad = arrgrad.view()
    part3grad.shape = (theta.size, gl.size, gl.size)

    for j in range(gl.numsteps-2):
        q += np.sum(part3 * allout[1][3][j])*gl.k*gl.k
        for i in range(n):
            gradq[i] = gradq[i] + np.sum(part3grad[i,...] * allout[1][3][j])*gl.k*gl.k

    return (-q),(-gradq)


def qfunopt(theta, allout, sp3):
    # compute part3 and its gradient which is then shared across processes
    mypairs = np.hstack((gl.pairs, np.repeat(a=np.array([theta]),repeats=gl.numpairs,axis=0)))
    arr = tonumpyarray(sp3[0])
    part3 = arr.view()
    part3.shape = (gl.size, gl.size)
    arrgrad = tonumpyarray(sp3[1])
    part3grad = arrgrad.view()
    part3grad.shape = (theta.size, gl.size, gl.size)
    logG(theta, part3, part3grad)

    with closing(mp.Pool(processes=1, initializer=initpart3, initargs=(sp3,))) as p:
        qfunall = p.map(qfun, zip(mypairs, allout))

    totq = 0.0
    totqgrad = np.zeros(len(theta))
#    for i in range(gl.numpairs):
#        totq += qfunall[i][0]
#        totqgrad += qfunall[i][1]
    for x in qfunall:
        totq += x[0]
        totqgrad += x[1]

    return totq, totqgrad


def Estep(shared_prop, thetak):
    arr = tonumpyarray(shared_prop)
    arrmat = arr.view()
    arrmat.shape = (gl.size, gl.size)
    integrandmat(thetak, arrmat)

    mypairs = np.hstack((gl.pairs, np.repeat(a=np.array([thetak]),repeats=gl.numpairs,axis=0)))

    with closing(mp.Pool(initializer=init, initargs=(shared_prop,))) as p:
        allout = p.map(dtqall, mypairs)

    # p.join()
    return list(allout)


def Mstep(sp3, thetak, allout):
    # goal is to minimize qfunopt(theta, allout) over theta
    # then we return the minimizer
    # qfunout = qfunopt(theta, allout)
    res = scipy.optimize.minimize(qfunopt, thetak, args=(allout, sp3), method='BFGS', jac=True, options={'disp':True})
    return res.x


def init(shared_prop_):
    global shared_prop
    shared_prop = shared_prop_


def initpart3(sp3_):
    global shared_part3
    global shared_part3grad
    shared_part3 = sp3_[0]
    shared_part3grad = sp3_[1]


def tonumpyarray(mp_arr):
    return np.frombuffer(mp_arr.get_obj())


def main():
    # logging
    # logger = mp.log_to_stderr()
    # logger.setLevel(logging.INFO)

    # create simulated data
    theta = np.array([1., 0.75, 2.])
    spdt = 1000
    idt = gl.deltat/spdt
    idt12 = np.sqrt(idt)
    sampsize = 24 + 1
    x = np.zeros(sampsize)
    t = np.zeros(sampsize)
    x[0] = 0.1
    t[0] = 0
    cur = x[0]
    curt = t[0]
    for i in range(sampsize-1):
        for j in range(spdt):
            cur += idt*f(theta,cur) + idt12*g(theta,cur)*np.asscalar(scipy.stats.norm.rvs(size=1))
            curt += idt
            
        x[i+1] = cur
        t[i+1] = curt

    # just for our purposes, print out max abs val
    print(np.max(np.abs(x)))

    # distribute the data in pairs plus counter
    gl.pairs = np.vstack((np.array([x[0:(len(x)-1)]]),np.array([x[1:len(x)]]))).T
    gl.numpairs = len(gl.pairs)

    thetak = np.array([0.5, 1.0, 1.5])

    # create shared memory array to store integrandmat
    shared_prop = mp.Array(ctypes.c_double, gl.size2)

    # store the grid matrices in shared memory
    gridxvec = tonumpyarray(gl.gridx)
    gridx = gridxvec.view()
    gridx.shape = (gl.size, gl.size)
    gridyvec = tonumpyarray(gl.gridy)
    gridy = gridyvec.view()
    gridy.shape = (gl.size, gl.size)
    gridx, gridy = np.meshgrid(gl.grid, gl.grid, sparse=False, indexing='ij')

    # create shared memory to store part3 logG matrix and its gradient
    shared_part3 = mp.Array(ctypes.c_double, gl.size2)
    shared_part3grad = mp.Array(ctypes.c_double, gl.size2 * thetak.size)
    sp3 = (shared_part3, shared_part3grad)

    # check gradient with finite difference
    # allout = Estep(thetak)
    # test = scipy.optimize.check_grad(objqfunopt, gradqfunopt, thetak, allout)
    # print(test)

    numiters = 100
    for i in range(numiters):
        print(i)
        print("E step")
        allout = Estep(shared_prop, thetak)
        print("M step")
        thetak[:] = Mstep(sp3, thetak, allout)
        print(thetak)

if __name__ == '__main__':
    main()

