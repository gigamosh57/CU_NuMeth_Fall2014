### Page Weil
### CVEN 5537
### 11/14/14
### Project 4, Problem 1

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

###################################
###### Begin define problem variables

tol = 10**-6

# Grid size
NX = 101
NY = 101
xmax = 1.
ymax = 1.
x = np.arange(0,xmax+xmax/(NX-1),xmax/(NX-1))
y = np.arange(0,ymax+ymax/(NY-1),ymax/(NY-1))

# initialize matrices
u = np.empty((NX,NY))
unew = u

# Boundary conditions
u[0,:]=x
u[NY-1,:]=x
u[:,0]=0
u[:,NX-1]=1

###### End define problem variables
###################################

###################################
###### Begin Jacobi solver



# while toleranceerror > tol
err = tol+1
while err > tol:
  # Loop through interior nodes
  for i in range(1,NX-1):
    for j in range(1,NY-1):
      delta = 1/4*(u[i-1,j]+u[i,j-1]+u[i,j+1]+u[i+1,j]-4*u[i,j])
      unew = u[i,j]+delta
  # calculate error
  err = np.max(delta)

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
CS = plt.contourf(X, Y, u, 10, # [-1, -0.1, 0, 0.1],
                        #alpha=0.5,
                        cmap=plt.cm.bone)#,
                        #origin=origin)

plt.title('Test')
plt.xlabel('x')
plt.ylabel('y')

#plt.figure()
plt.show()

###### End plotting code
###################################


