### Page Weil
### CVEN 5537
### 11/14/14
### Project 4, Problem 1c - Vectorized Jacobi

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

###################################
###### Begin define problem variables

tol = 10**-6
maxiter = 500000

# Grid size
NX = 101
NY = 101

xmax = 1.
ymax = 1.
x = np.arange(0,xmax+xmax/(NX-1),xmax/(NX-1))
y = np.arange(0,ymax+ymax/(NY-1),ymax/(NY-1))

# initialize matrices
u = np.empty((NX,NY))*0.
unew = np.empty((NX,NY))*0.
delta = np.empty((NX,NY))*0.

# Boundary conditions
u[0,:]=x
u[NY-1,:]=x
u[:,0]=0.
u[:,NX-1]=1.

###### End define problem variables
###################################

start1 = datetime.now()
###################################
###### Begin Jacobi solver

# while toleranceerror > tol
err = tol+1
it=0
while err > tol:
  it = it + 1
  
  # define lookup vectors
  i = np.arange(1,NX-1)
  j = np.arange(1,NY-1)
  # vectorize delta
  
  ### THIS DOESN'T WORK YET
  # define range based on matrix layout
  # convert delta back to right array shape once calculated
  delta = 1/4.*(u[i-1,j]+u[i,j-1]+u[i,j+1]+u[i+1,j]-4.*u[i,j])
  ### THIS DOESN'T WORK YET
  
  unew =  u+delta
  u = unew
  # store in u
  
  # calculate error
  err = np.max(delta)
  if it > maxiter: 
    err = tol
    print("Hit ",str(maxiter)," iter")


# Set special conditions (Type II BCs?)
# Set values for exterior nodes

###### End Jacobi solver
###################################

end1 = datetime.now()

print("Vectorized Jacobi took ",str(it)," iterations and "+str(end1-start1))


