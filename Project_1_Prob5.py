#### Page Weil
#### 9/15/14
#### CVEN 5537
#### PROJECT 1, PROBLEM 5

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pylab as plt

thetasum = lambda x, n, t : (np.sin(((2*n+1)*np.pi*x)/(2*L))*np.exp(-1*(((2*n+1)**2)*(np.pi**2)*alpha*t)/(4*L**2)))/((2*n+1)*math.pi)

#### Part (a), use for loops

## Initialize default values
L = 1
alpha = 1
nmax = 50
theta0 = 1

## Initialize arrays
t = [0.01,0.1,0.2,0.5,1,2,5]
n = np.arange(0,nmax+1)
dx = 0.05
xeval = list(np.arange(0,L+dx,dx))
thetas = np.zeros((len(t),len(xeval)))

## Loop through and solve
for tn in t:
	tx = t.index(tn)
	for xn in xeval:
		xx = xeval.index(xn)
		thetanow = 0
		for nn in n:
			thetanow = thetanow + thetasum(xn,nn,tn)
		thetas[tx,xx] = theta0*(1-4*thetanow)

## Loop through and plot
for p in np.arange(0,len(thetas[:,1])):
	plt.plot(xeval,thetas[p,:])

plt.show()


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

