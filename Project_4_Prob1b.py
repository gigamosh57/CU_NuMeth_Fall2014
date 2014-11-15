### Page Weil
### CVEN 5537
### 11/14/14
### Project 4, Problem 1b - GS/SOR

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

k=2.975
w = 0.5*(1-k/NX)

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

###################################
###### Begin GS/SOR solver

# while toleranceerror > tol
err = tol+1
it=0
while err > tol:
  it = it + 1
  #print(str(iter))
  # Loop through interior nodes
  for i in range(1,NX-1):
    for j in range(1,NY-1):
        delta[i,j] = w*(u[i-1,j]+u[i,j-1]+u[i,j+1]+u[i+1,j]-4.*u[i,j])
        u[i,j] = u[i,j]+delta[i,j]
  # calculate error
  err = np.max(delta)
  if it > maxiter: 
    err = tol
    print("Hit ",str(maxiter)," iter")

print("For k = ",str(k)," SOR took ",str(it)," iterations")

###### End GS/SOR solver
###################################

