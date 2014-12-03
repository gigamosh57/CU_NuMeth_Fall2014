### Page Weil
### CVEN 5537
### 11/14/14
### Project 4, Problem 1a - Jacobi

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
u = np.zeros((NX,NY))
unew = np.zeros((NX,NY))
delta = np.zeros((NX,NY))

# Boundary conditions
u[0,:]=x
u[NY-1,:]=x
u[:,0]=0.
u[:,NX-1]=1.

###### End define problem variables
###################################

###################################
###### Begin Jacobi solver

# while toleranceerror > tol
err = tol+1
it=0
while err > tol:
  it = it + 1
  #print(str(iter))
  # Loop through interior nodes
  for i in range(1,NX-1):
    for j in range(1,NY-1):
      delta[i,j] = 1/4.*(u[i-1,j]+u[i,j-1]+u[i,j+1]+u[i+1,j]-4.*u[i,j])
#      unew[i,j] = u[i,j]+delta[i,j]
  unew =  u+delta
  u = unew
  # calculate error
  err = np.max(delta)
  if it > maxiter: 
    err = tol
    print("Hit ",str(maxiter)," iter")

print("Jacobi took ",str(it)," iterations")
# Set special conditions (Type II BCs?)
# Set values for exterior nodes


###### End Jacobi solver
###################################


###################################
###### Begin plotting code here

# Make plottable arrays (X,Y)
# Coordinate grids 
X, Y = np.meshgrid(x, y)

# use this to pull in u values into plottable matrix
#   itemindex = numpy.where(array==item)

# We are using automatic selection of contour levels;
# this is usually not such a good idea, because they don't
# occur on nice boundaries, but we do it here for purposes
# of illustration.
CS = plt.contourf(X, Y, u, 10, cmap=plt.cm.bone)

plt.title('Jacobi')
plt.xlabel('x')
plt.ylabel('y')

#plt.figure()
plt.show()

###### End plotting code
###################################


