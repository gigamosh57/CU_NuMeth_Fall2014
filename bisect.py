## This code calculates the roots of an equation from prob 5, hw 2
import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt
# for fsolve
from scipy.optimize import fsolve

# Function
f = lambda x : x * np.tan(x) - r

for r in np.arange(0,10.5,0.5):
  print proot(bd
# Bisection solver begins here
bd = [0,1.57]
r = 3
def proot(bounds,n = 1000,tol = 0.01): 
  it = 0
  i = 0
  while (i < n):
      i = i + 1
      it = i
      xmid = np.mean(bounds)
      t = [f(xmid),f(bounds[0])] 
      if (all( item > 0 for item in t) or all( item <= 0 for item in t)):
            bounds[0] = xmid
      else: 
           t = [f(xmid),f(bounds[1])] 
           if (all( item > 0 for item in t) or all( item <= 0 for item in t)):
                bounds[1] = xmid
      if (np.absolute(f(xmid)) < tol):
            i = n+1
  return it, xmid, f(xmid)

### LOOPED SOLUTION HERE


iter = ()
tanroot = ()
root = ()
rplot = ()
for r in np.arange(0,10.5,0.5):
  bd = [0,1.57]
  a,b,c = proot(bd)
  rplot = rplot + (r,)
  iter = iter + (a,)
  tanroot = tanroot + (b,)
  root = root + (c,)

print rplot
print tanroot
print iter

plt.plot(rplot,iter)
plt.plot(rplot,tanroot)
plt.show()
