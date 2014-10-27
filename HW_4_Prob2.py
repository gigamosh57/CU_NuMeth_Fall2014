#### Page Weil
#### 10/27/14
#### CVEN 5537
#### HOMEWORK 4, PROBLEM 2

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pylab as plt



# OLD, REMOVE  thetasum = lambda x, n, t : (np.sin(((2*n+1)*np.pi*x)/(2*L))*np.exp(-1*(((2*n+1)**2)*(np.pi**2)*alpha*t)/(4*L**2)))/((2*n+1)*math.pi)

#### Part (a), Example 1, analytical solution to heat equation for specific BCs

usum = lambda x, m, t: 2*u0/m*np.pi*(1*(-1)**m)*np.sin(m*np.pi*x/L)*np.exp(-1*m**2*np.pi**2*b*t/L**2)

## Initialize default values
L = 1
b = 1
mmax = 50
u0 = 1

## Initialize arrays
t = [0.01,0.1,0.2,0.5,1,2,5]
m = np.arange(0,nmax+1)
dx = 0.05
xeval = list(np.arange(0,L+dx,dx))
u_values = np.zeros((len(t),len(xeval)))

## Loop through and solve
for tm in t:
	tx = t.index(tm)
	for xm in xeval:
		xx = xeval.index(xm)
		unow = 0
		for mm in m:
			unow = unow + usum(xm,mm,tm)
		u_values[tx,xx] = unow

## Loop through and plot
fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
axes.set_xlabel('x (L)')
axes.set_ylabel('u')
for p in np.arange(0,len(u_values[:,1])):
	axes.plot(xeval,u_values[p,:],color = colors[p],label = "t: " + str(t[p]))

axes.legend()
plt.ylim((0,1))
plt.show()


plt.show()

#### Part (b), Example 2, analytical solution to heat equation for specific BCs

#### First find the appropriate am values for this solution
