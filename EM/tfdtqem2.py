
# coding: utf-8

# In[13]:

import tensorflow as tf


# In[27]:

import numpy as np


# In[50]:

class gl:
    myftype = 'float64'
    myitype = 'int64'
    k = tf.constant(0.05,dtype=myftype)
    M = tf.constant(80,dtype=myitype)
    deltat = tf.constant(0.1,dtype=myftype)
    numsteps = 3
    h = tf.divide(deltat,numsteps)
    grid = tf.multiply(k,tf.cast(tf.range(tf.negative(M), tf.add(M, 1), delta=1),dtype=myftype))
    gridx, gridy = tf.meshgrid(grid, grid, indexing='ij')
    gridsize = tf.add(tf.multiply(M, 2), 1)
    gridsize2 = tf.pow(gridsize, 2)


# In[51]:

def tff(theta, x):
    fval = tf.multiply(theta[0],tf.subtract(theta[1],x))
    return fval

def tfg(theta, x):
    gval = tf.multiply(theta[2],tf.ones(tf.shape(x),dtype=gl.myftype))
    return gval


# In[52]:

def integrandmat(inx, iny, th):
    my2 = tf.constant(2.0,gl.myftype)
    tfmu = tf.add(iny,tf.multiply(tff(theta=th,x=iny),gl.h))
    tfsig = tf.multiply(tf.sqrt(gl.h),tfg(theta=th,x=iny))
    tfc0 = tf.reciprocal(tf.multiply(tf.sqrt(tf.multiply(my2,tf.constant(np.pi,dtype=gl.myftype))),tfsig))
    tfnumer = tf.negative(tf.square(tf.subtract(inx,tfmu)))
    tfdenom = tf.multiply(my2,tf.square(tfsig))
    tfprop = tf.multiply(tfc0,tf.exp(tf.divide(tfnumer,tfdenom)))
    return tfprop


# In[53]:

def dtqall(theta, init, final):
    # supposed to be a column vector
    lamb = tf.expand_dims(integrandmat(gl.grid, init, theta),1)
    
    # supposed to be a row vector
    gamm = gl.k * tf.expand_dims(integrandmat(final, gl.grid, theta),0)
    
    # declare two tensors in advance
    
    lamblist = []
    gammlist = []
    lamblist.append(lamb)
    gammlist.append(gamm)
    
    for j in range(gl.numsteps-2):
        lamblist.append(gl.k*tf.matmul(gl.A,lamblist[j]))
        gammlist.append(gl.k*tf.matmul(gammlist[j],gl.A))
        
    complete = tf.reshape(tf.matmul(gammlist[0],lamblist[gl.numsteps-2]),shape=[])
    
    # both first and last wil be column vectors
    c = tf.multiply(gl.k,complete)
    first = tf.divide(tf.multiply(tf.transpose(gammlist[gl.numsteps-2]),lamblist[0]),c)
    last = tf.divide(tf.multiply(tf.transpose(gammlist[0]),lamblist[gl.numsteps-2]),c)

    # compute all intermediate 2d pdfs
    del lamblist[-1]
    del gammlist[-1]
    lambtensor = tf.stack(lamblist)
    # diag = []
    # diag.append(tf.shape(lambtensor))
    gammtensor = tf.stack(gammlist[::-1])
    # diag.append(tf.shape(gammtensor))
    A0 = tf.expand_dims(gl.A,0)
    # diag.append(tf.shape(A0))
    pdf2dlist = tf.transpose(tf.multiply(gammtensor, lambtensor),perm=[0,2,1])
    # diag.append(tf.shape(pdf2dlist))
    pdf2dlist = tf.multiply(pdf2dlist, A0)
    pdf2dlist = tf.divide(pdf2dlist, c)

    return complete, first, last, pdf2dlist #, diag


# In[54]:

def logG(x, y, theta):
    fv = tff(theta,y)
    gv = tfg(theta,y)
    mu = tf.add(y,tf.multiply(fv,gl.h))
    pr = tf.subtract(x,mu)
    pr2 = tf.square(pr)
    gv2 = tf.square(gv)
    my2 = tf.constant(2.0,dtype=gl.myftype)
    mypi = tf.constant(np.pi,dtype=gl.myftype)
    lgp1 = tf.negative(tf.divide(tf.log(tf.multiply(my2*mypi*gl.h,gv2)),my2))
    lgp2 = tf.negative(tf.divide(pr2,tf.multiply(my2*gl.h,gv2)))
    lg = tf.add(lgp1,lgp2)        
    return lg


# In[55]:

# allout is the result of a call to dtqall
# here x is the data
def qfun(theta, allout, init, final):
    q = []

    # first term in the summation (j=1 case)
    part1 = logG(gl.grid,init,theta)
    q.append(tf.reduce_sum(tf.multiply(part1,allout[1][:,0]))*gl.k)

    # last term in the summation (j=F case)
    part2 = logG(final,gl.grid,theta)
    q.append(tf.reduce_sum(tf.multiply(part2,allout[2][:,0]))*gl.k)
    
    # all intermediate terms
    part3 = logG(gl.gridx,gl.gridy,theta)
    # for j in range(gl.numsteps-2):
    #     q.append(tf.tensordot(part3,allout[3][j,:],axes=[[0,1],[0,1]])*gl.k*gl.k)
    # test = tf.add_n(test)
    
    q.append(tf.reduce_sum(tf.multiply(tf.expand_dims(part3,0),allout[3]))*gl.k*gl.k)
    
    qout = tf.negative(tf.add_n(q))

    return qout


# In[ ]:



# print(sess.run(qfun(th, allout, init, final), feed_dict = {th : thval}))
# print(sess.run(tf.gradients(qfun(th, allout, init, final),th), feed_dict = {th : thval}))
# part3 = logG(gl.gridx,gl.gridy,th)
# print(sess.run(tf.shape(allout[3][0,:]), feed_dict = {th : thval}))
# print(sess.run(tf.tensordot(part3,allout[3][0,:],axes=[[0,1],[0,1]]), feed_dict = {th : thval}))
# print(sess.run(fir, feed_dict = {th : thval}))
# print(sess.run(gl.k*gl.k*tf.reduce_sum(pdf,axis=[1,2]), feed_dict = {th : thval}))
# print(sess.run(tf.gradients(logG(gl.gridx,gl.gridy,th),th),  feed_dict = {th : thval}))


# In[78]:

# load simulated data from the file
import csv
datain = []
with open('simdata.csv', newline='') as f:
    reader = csv.reader(f,quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        datain.append(row)
        
xt = np.array(datain,dtype=gl.myftype)


# In[113]:

nrow = xt.shape[0]
pairs = np.vstack([xt[0:(nrow-1),0],xt[1:nrow,0]]).T


# In[ ]:

# start with some initial theta, call that thetak
thval = [0.5,0.7,0.2]
def dtqallwrap(onepair):
    return dtqall(theta=gl.thetak,init=onepair[0],final=onepair[1])

def myloss(th, allallouts, pairs):
    newlist = []
    for i in range(len(pairs)):
        newlist.append([allallouts[i],pairs[i]])
    
    temp = [qfun(th, x[0], x[1][0], x[1][1]) for x in newlist]
    qout = tf.reduce_sum(temp)
    # qout = tf.Variable(tf.constant(0.,dtype=gl.myftype))
    #qout = tf.constant(0.,dtype=gl.myftype)
    #for i in range(len(pairs)):
    #    qout += qfun(th,allallouts[i],pairs[i][0],pairs[i][1])
        
    return qout

emiters = 100
for emi in range(emiters):
    sess = tf.Session()

    # E step
    gl.thetak = tf.constant(thval,dtype=gl.myftype)
    gl.A = integrandmat(gl.gridx, gl.gridy, gl.thetak)
    allallouts = sess.run(list(map(dtqallwrap, pairs)))

    # M step
    # define the actual q function
    th = tf.Variable(gl.thetak,dtype=gl.myftype)
    sess.run(tf.global_variables_initializer())
    loss = myloss(th, allallouts, pairs)
    print("Loss before M step: " + str(sess.run(loss)))

    # optimize the loss
    if emi==0:
        learning_rate = 0.001
        numiters = 1
    if emi>=1:
        learning_rate = 0.01
        numiters = 100

    my_opt = tf.train.GradientDescentOptimizer(learning_rate)
    train_step = my_opt.minimize(loss)
    for i in range(numiters):
        sess.run(train_step)

    print("Loss after M step: " + str(sess.run(loss)))

    thval = sess.run(th)
    sess.close()


