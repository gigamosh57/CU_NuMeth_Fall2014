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

# Grid definition
#deltax = 0.01
nx = int(L/deltax+1)
deltat = 0.5*deltax**2
tmax = 0.5
mu = b*deltat/deltax**2

# Initialize variables
a = np.arange(0,nx)*0.0
b = np.arange(0,nx)*0.0
c = np.arange(0,nx)*0.0
r = np.arange(0,nx)*0.0

#BC Definition
c1 = 1.0 # T at x = 0
c2 = 0.0 # dT/dx at x = 0
c3 = 0.0 # T at x = L
c4 = 1.0 # dT/dx at x = L

# IC definition
x = np.arange(0,nx)*0.0
x[0] = 1

# for each TS, determine a, b, c and r

#a
a[0] = 0
a[1:(nx-1)] = -mu/2
a[nx-1] = -c4/deltax

# b
b[0] = c1 - c2/deltax
b[1:(nx-1)] = 1+mu
b[nx-1] = c3 + c4/deltax
 
# c
c[0] = c2/deltax
c[1:(nx-1)] = -mu/2
c[nx-1] = 0
for t in np.arange(0,tmax,deltat):
   
   # for each TS, determine a, b, c and r
   
   # r
   r[0] = c1*x[0] + c2*(x[1]-x[0])/(deltax)
   r[1:(nx-1)] = (mu/2)*x[0:(nx-2)] + (1-mu)*x[1:(nx-1)] + (mu/2)*x[2:nx]
   r[nx-1] = c3*x[0] + c4*(x[J-1]-x[J-2])/(deltax)
   
   ####### THOMAS ALGORITHM HERE
   # plug into Thomas
   x = thomas(a,b,c,r)

xplot = np.arange(0,L+deltax,deltax)
axes.plot(xplot,x,color = [rand.random(),rand.random(),rand.random()],label = "dt: " + str(deltat))
case1 = x[nx/2]

######## End Problem Part a
################################################

################################################
######## Setup Problem Part b

# Variable definition
L = 1.0
b = 1.0

# Grid definition
#deltax = 0.01
nx = int(L/deltax+1)
deltat = 5*deltax**2
#tmax = 0.5
mu = b*deltat/deltax**2

# Initialize variables
a = np.arange(0,nx)*0.0
b = np.arange(0,nx)*0.0
c = np.arange(0,nx)*0.0
r = np.arange(0,nx)*0.0

#BC Definition
c1 = 1.0 # T at x = 0
c2 = 0.0 # dT/dx at x = 0
c3 = 0.0 # T at x = L
c4 = 1.0 # dT/dx at x = L
f1 = 1.0 # 
f2 = 1.0 #

# IC definition
x = np.arange(0,nx)*0.0
x[0] = 1

#a
a[0] = 0
a[1:(nx-1)] = -mu/2
a[nx-1] = -c4/deltax

# b
b[0] = c1 - c2/deltax
b[1:(nx-1)] = 1+mu
b[nx-1] = c3 + c4/deltax
 
# c
c[0] = c2/deltax
c[1:(nx-1)] = -mu/2
c[nx-1] = 0
for t in np.arange(0,tmax,deltat):
   
   # for each TS, determine a, b, c and r
   
   # r
   r[0] = c1*x[0] + c2*(x[1]-x[0])/(deltax)
   r[1:(nx-1)] = (mu/2)*x[0:(nx-2)] + (1-mu)*x[1:(nx-1)] + (mu/2)*x[2:nx]
   r[nx-1] = c3*x[0] + c4*(x[J-1]-x[J-2])/(deltax)
   
   ####### THOMAS ALGORITHM HERE
   # plug into Thomas
   x = thomas(a,b,c,r)

xplot = np.arange(0,L+deltax,deltax)
axes.plot(xplot,x,color = [rand.random(),rand.random(),rand.random()],label = "dt: " + str(deltat))
case2 = x[nx/2]
######## End Problem Part b
################################################

################################################
######## Setup Problem Part c

# Variable definition
L = 1.0
b = 1.0

# Grid definition
#deltax = 0.01
nx = int(L/deltax+1)
deltat = 50*deltax**2
#tmax = 0.5
mu = b*deltat/deltax**2

# Initialize variables
a = np.arange(0,nx)*0.0
b = np.arange(0,nx)*0.0
c = np.arange(0,nx)*0.0
r = np.arange(0,nx)*0.0

#BC Definition
c1 = 1.0 # T at x = 0
c2 = 0.0 # dT/dx at x = 0
c3 = 0.0 # T at x = L
c4 = 1.0 # dT/dx at x = L

# IC definition
x = np.arange(0,nx)*0.0
x[0] = 1

#a
a[0] = 0
a[1:(nx-1)] = -mu/2
a[nx-1] = -c4/deltax

# b
b[0] = c1 - c2/deltax
b[1:(nx-1)] = 1+mu
b[nx-1] = c3 + c4/deltax
 
# c
c[0] = c2/deltax
c[1:(nx-1)] = -mu/2
c[nx-1] = 0
for t in np.arange(0,tmax,deltat):
   
   # for each TS, determine a, b, c and r
   
   # r
   r[0] = c1*x[0] + c2*(x[1]-x[0])/(deltax)
   r[1:(nx-1)] = (mu/2)*x[0:(nx-2)] + (1-mu)*x[1:(nx-1)] + (mu/2)*x[2:nx]
   r[nx-1] = c3*x[0] + c4*(x[J-1]-x[J-2])/(deltax)
   
   ####### THOMAS ALGORITHM HERE
   # plug into Thomas
   x = thomas(a,b,c,r)

xplot = np.arange(0,L+deltax,deltax)
axes.plot(xplot,x,color = [rand.random(),rand.random(),rand.random()],label = "dt: " + str(deltat))
case3 = x[nx/2]
######## End Problem Part c
################################################

################################################
######## Begin Analytical Solution

#### Page Weil
#### 9/15/14
#### CVEN 5537
#### PROJECT 1, PROBLEM 5

thetasum = lambda x, n, t : (np.sin(((2*n+1)*np.pi*x)/(2*L))*np.exp(-1*(((2*n+1)**2)*(np.pi**2)*alpha*t)/(4*L**2)))/((2*n+1)*math.pi)

#### Part (b), do summation with matrix capabilities

## Initialize default values
L = 1
alpha = 1
nmax = 100
theta0 = 1

## Initialize arrays
t = [0.01,tmax]

n = np.arange(0,nmax+1)
dx = deltax
xeval = list(np.arange(0,L+dx,dx))
thetas = np.zeros((len(t),len(xeval)))

## Loop through and solve
for tn in t:
   tx = t.index(tn)
   for xn in xeval:
      xx = xeval.index(xn)
      thetas[tx,xx] = theta0*(1-4*np.sum(thetasum(xn,n,tn)))

axes.plot(xeval,thetas[1,:],color = [rand.random(),rand.random(),rand.random()],label = "ana, t:"+str(t[1]))
analytical = thetas[1,int(nmax/2)]

axes.legend()
plt.show()

print("For x=0.5, t=0.5; Case 1, x="+
   str(round(case1,8))+" Case 2, x="+
   str(round(case2,8))+" Case 3, x="+
   str(round(case3,8))+" Analytical, x="+
   str(round(analytical,8)))

