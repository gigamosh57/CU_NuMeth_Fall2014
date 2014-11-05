### Page Weil
### CVEN 5537
### 11/3/14
### Project 3, Problem 2

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pylab as plt

################################################
######## Thomas Algorithm
def thomas(a,b,c,r): 
   # Define vectors as size of b
   e = b*0
   f = b*0
   g = b*0
   y = b*0
   x = b*0
   J = len(b)
  
   # define efg matrix
   f[0] = b[0]
   g[0] = c[0]/f[0]
   ########################################
   ######## PROBLEM IS HERE IN HOW F IS DEFINED
   ########################################
   for i in np.arange(1,J):
      e[i] = a[i]
      f[i] = b[i]-e[i]*g[i-1]
      g[i] = c[i]/f[i]
   
   # solve for y
   y[0] = r[0]/f[0] 
   for i in np.arange(1,J-1):
      y[i] = (r[i]-e[i]*y[i-1])/f[i]
   
   # solve for x
   x[J-1] = y[J-1]
   for i in np.arange(0,J-2)[::-1]:
      x[i] = (y[i]-g[i]*x[i+1])
      
   return x
################################################
######## Crank-Nicholson implementation

# Problem definition
L = 1
b = 1

#BC Definition
c1 = 1 # T at x = 0
c2 = 0 # dT/dx at x = 0
c3 = 0 # T at x = L
c4 = 1 # dT/dx at x = L

# Grid definition
deltax = 0.01
nx = int(L/deltax+1)
deltat = 0.5*np.sqrt(deltax)
tmax = 0.5
mu = b*deltat/deltax**2

a = np.arange(0,nx)
b = np.arange(0,nx)
c = np.arange(0,nx)
r = np.arange(0,nx)

# IC definition
x = np.arange(0,nx)
x[0] = 1

for t in np.arange(0,tmax,deltat):
   # for each TS, determine a, b, c and r
   
   # a
   a[0] = 0
   a[1:(nx-1)] = -mu
   a[nx-1] = -c4/deltax
   
   # b
   b[0] = c1 - c2/deltax
   b[1:(nx-1)] = 1+2*mu
   b[nx-1] = c3-c4/deltax
   
   # c
   c[0] = c2/deltax
   c[1:(nx-1)] = -mu
   c[nx-1] = 0
   
   # plug into Thomas
   x = thomas(a,b,c,r)
   
   # r
   r[0] = c1*x[0] + c2*(x[1]-x[0])/deltax
   r[1:(nx-2)] = x[1:(nx-2)]
   r[nx-1] = c3*x[0] + c4*(x[J-1]-x[J-2])/deltax

################################################
######## Project 1, Problem 5 solution for comparison

thetasum = lambda x, n, t : (np.sin(((2*n+1)*np.pi*x)/(2*L))*np.exp(-1*(((2*n+1)**2)*(np.pi**2)*alpha*t)/(4*L**2)))/((2*n+1)*math.pi)

#### Part (b), do summation with matrix capabilities

## Initialize default values
L = 1
alpha = 1
nmax = 50
theta0 = 1

## Initialize arrays
t = [0.01,0.1,0.2,0.5,1,2,5]
colors = ()
for a in np.arange(0,len(t)): colors = colors + ([(float(a)/(len(t))),0,0],)

n = np.arange(0,nmax+1)
dx = 0.05
xeval = list(np.arange(0,L+dx,dx))
thetas = np.zeros((len(t),len(xeval)))

## Loop through and solve
for tn in t:
   tx = t.index(tn)
   for xn in xeval:
      xx = xeval.index(xn)
      thetas[tx,xx] = theta0*(1-4*np.sum(thetasum(xn,n,tn)))

## Loop through and plot
fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
axes.set_xlabel('x (L)')
axes.set_ylabel('theta')
for p in np.arange(0,len(thetas[:,1])):
   axes.plot(xeval,thetas[p,:],color = colors[p],label = "t: " + str(t[p]))

axes.legend()
plt.ylim((0,1))
plt.show()

### Both plots show the same result.  Since x = 0 is kept at a constant temp, it makes sense
### that all the lines converge on x = 0, theta = 1.  As time goes on, the whole rod slowly heats
### with the portion closer to the heat source heating faster than the portion that is not insulated

