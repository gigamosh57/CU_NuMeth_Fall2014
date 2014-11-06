#thomas

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pylab as plt
# for colors
import random as rand

################################################
######## Thomas Solver

def thomas(a,b,c,r):
   # Define vectors as size of b
   e = b*0
   f = b*0
   g = b*0
   y = b*0
   J = len(b)
   
   # r
   r[0] = c1*x[0] + c2*(x[1]-x[0])/deltax
   r[1:(nx-1)] = x[1:(nx-1)]
   r[nx-1] = c3*x[0] + c4*(x[J-1]-x[J-2])/deltax
   
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

################################################
######## Initialize figure
fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
axes.set_xlabel('x (L)')
axes.set_ylabel('T')
######## End Initialize figure
################################################


################################################
######## Setup Problem Part a

# Variable definition
L = 1.0
b = 1.0
h0 = 2
t0 = 1

# Grid definition
deltax = 0.01
nx = int(L/deltax+1)
deltat = 0.5*np.sqrt(deltax)
tmax = 0.5
mu = b*deltat/deltax**2

# Initialize variables
a = np.arange(0,nx)*0.0
b = np.arange(0,nx)*0.0
c = np.arange(0,nx)*0.0
r = np.arange(0,nx)*0.0

#BC Definition
c1 = 0.0 # T at x = 0
c2 = 0.0 # dT/dx at x = 0
c3 = 0.0 # T at x = L
c4 = 0.0 # dT/dx at x = L
f1 = 1.0 # 
f2 = 1.0 #

# IC definition
x = np.arange(0,nx)*0.0
x[0] = h0


for t in np.arange(0,tmax,deltat):
  # Changing boundary conditions from problem definition
  c1 = h0*np.sin(2*np.pi*t/t0)
  
  # Determine a, b, c and r at the outset
  #a
  a[0] = 0
  a[1:(nx-1)] = -mu
  a[nx-1] = -c4/deltax
  
  # b
  b[0] = c1 - c2/deltax
  b[1:(nx-1)] = 1+2*mu
  b[nx-1] = c3 + c4/deltax
  
  # c
  c[0] = c2/deltax
  c[1:(nx-1)] = -mu
  c[nx-1] = 0
  ####### THOMAS ALGORITHM HERE
  # plug into Thomas
  x = thomas(a,b,c,r)

xplot = np.arange(0,L+deltax,deltax)
axes.plot(xplot,x,color = [rand.random(),rand.random(),rand.random()],label = "dt: " + str(deltat))

######## End Problem Part a
################################################

axes.legend()
plt.show()


