### Page Weil
### CVEN 5537
### 11/14/14
### Final Project, Prob 1

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pylab as plt
# for colors
import random as rand
# for timing
from datetime import datetime

################################################
######## Thomas Solver

def thomas(a,b,c,r):
   # Define vectors as size of b
   e = b*0
   f = b*0
   g = b*0
   y = b*0
   J = len(b)
   
   # define efg matrix
   f[0] = b[0]
   g[0] = c[0]/f[0]
   for i in np.arange(1,J):
      e[i] = a[i]
      f[i] = b[i]-e[i]*g[i-1]
      g[i] = c[i]/f[i]
   
   # solve for y
   y[0] = r[0]/f[0] 
   for i in np.arange(1,J):
      y[i] = (r[i]-e[i]*y[i-1])/f[i]
   
   # solve for x
   x[J-1] = y[J-1] 
   for i in np.arange(0,J-1)[::-1]:
      x[i] = (y[i]-g[i]*x[i+1])
   return x

######## End Thomas Solver
################################################


###################################
###### Begin define problem variables

tol = 10**-6
maxiter = 100

# Grid size
dx = 0.01
dt = 0.01
xmax = 1
tmax = 100

NX = int(xmax/dx)+1
NT = int(tmax/dt)+1

x = np.arange(0,xmax+dx,dx)
t = np.arange(0,tmax+dt,dt)

alpha = 1.
p = 0 # 0,1,2,4

mu = alpha*dt/dx**2

# Boundary conditions
c1 = 1
c2 = 0
c3 = 1
c4 = 0
b1 = 0
bj = 1

# initialize matrices and initial conditions 
u = np.zeros((NX,NT)) # starts at zero everywhere except the end
u[NX-1,:]=bj
unew = np.zeros((NX)) 
uold = np.zeros((NX))  
a = np.zeros((NX)) 
a[NX-1]=-c4/dx
b = np.zeros((NX)) 
b[0] = c1-c2/dx
b[NX-1] = c3+c4/dx
c = np.zeros((NX)) 
c[0] = c2/dx
f = np.zeros((NX)) 

###### End define problem variables
###################################

# time loop 
for ti in range(1,NT):
   unew = u[:,ti-1]
   uold = u[:,ti-1]
   # loop through evaluations of [F]
   err = errtol+1
   while err > errtol:
      # check convergence
      err = np.amax(np.absolute(f))
      # evaluate f
      f[0] = (c1-c2/dx)*unew[0]+(c2/dx*unew[1])-b1
      f[1:(NX-1)] = -mu/(p+1)*(unew[0:NX-2]**(p+1))+(unew[1:NX-1]+2*mu/(p+1)*unew[1:NX-1]**(p+1))-mu/(p+1)*(unew[2:NX]**(p+1))-uold[1:NX-1]
      f[NX-1] = (-c4/dx)*unew[NX-2]+(c3+c4/dx)*unew[NX-1]-bj
      # create jacobian (as abc)
      a[1:NX] = -mu*(unew[0:(NX-1)])**p
      b[1:(NX-1)] = 1+2*mu*(unew[1:(NX-1)])**p
      c[0:(NX-1)] = -mu*(unew[1:(NX)])**p
      # solve with thomas
      delu = thomas(a,b,c,-1*f)
      # update value
      unew = unew + delu
      print(str(err))
    
  # Update u with converged values
  u[[:,ti]=unew
  print("IT: ",str(ti))
  
# plot multiple time values on single plot for single value of p
