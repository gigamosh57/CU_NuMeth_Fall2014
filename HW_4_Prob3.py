#### Page Weil
#### 10/28/14
#### CVEN 5537
#### HOMEWORK 4, PROBLEM 3

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pylab as plt

am = lambda m: -4*np.cos(np.pi*m)/(np.pi*(2*m+1))
alpha = lambda m: ((2*m+1)*np.pi/(2*L))
vsum = lambda x, m, t: am(m)*np.cos(alpha(m)*x)*np.exp(-1*(alpha(m)**2)*b*t)
uss = lambda x: 1
## Initialize default values
L = 1
b = 1
mmax = 50

## Initialize arrays
t = [0.01,0.05,0.2,0.3,0.5,1,2]
m = np.arange(0,mmax+1)
dx = 0.02
xeval = list(np.arange(0,L+dx,dx))
u_values = np.zeros((len(t),len(xeval)))

## Loop through and solve
for tm in t:
 tx = t.index(tm)
 for xm in xeval:
  xx = xeval.index(xm)
  vnow = 0
  for mm in m:
   vnow = vnow + vsum(xm,mm,tm)
  u_values[tx,xx] = vnow + uss(xm)

## Loop through and plot
colors = ()
for a in np.arange(0,len(t)): colors = colors + ([(float(a)/(len(t))),0,0],)

fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
axes.set_xlabel('x (L)')
axes.set_ylabel('u')
for p in np.arange(0,len(u_values[:,1])):
 axes.plot(xeval,u_values[p,:],color = colors[p],label = "t: " + str(t[p]))

axes.legend()
#plt.ylim((0,1))

plt.show()
#### First find the appropriate am values for this solution
