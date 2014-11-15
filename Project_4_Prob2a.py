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
maxiter = 5000
w = 1/4.

# Grid size
NX = 21
NY = 21

xmax = 1.
ymax = 1.
x = np.arange(0,xmax+xmax/(NX-1),xmax/(NX-1))
y = np.arange(0,ymax+ymax/(NY-1),ymax/(NY-1))

# initialize matrices
u = np.zeros((NX,NY))
unew = np.zeros((NX,NY))
delta = np.zeros((NX,NY))

# Boundary conditions
u[0,:]=0.
u[NY-1,:]=0.
u[:,0]=0.
u[:,NX-1]=0.

# set of roof points

x1 = int(NX/4)
x2 = int(NX*3/4)
y1 = int(NY/4)
y2 = int(NY/2)
y3 = int(NY*3/4)

pts = np.array([[x1,y1],[x1,y2],[x1,y3],[x2,y1],[x2,y2],[x2,y3]])
for a in pts:
   u[a[0],a[1]] = 1
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
	  if [i,j] not in pts.tolist():
		delta[i,j] = w*(u[i-1,j]+u[i,j-1]+u[i,j+1]+u[i+1,j]-4.*u[i,j])
	  u[i,j] = u[i,j]+delta[i,j]
  # calculate error
  err = np.max(delta)
  if it > maxiter: 
    err = tol
    print("Hit ",str(maxiter)," iter")

print("GS took ",str(it)," iterations")
# Set special conditions (Type II BCs?)
# Set values for exterior nodes


###### End GS/SOR solver
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
CS = plt.contourf(X, Y, u, 20, cmap=plt.cm.bone)

plt.title('Test')
plt.xlabel('x')
plt.ylabel('y')

#plt.figure()
plt.show()

###### End plotting code
###################################


