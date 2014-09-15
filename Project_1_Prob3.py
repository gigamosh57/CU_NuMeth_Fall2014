#### NUMERICAL METHODS PROJECT 1, PROBLEM 3
import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt

######## Bisectionsolver begins here
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
######## Newton solver ends here



# Define Function
f = lambda x : x*np.tan(x)-3
# Fprime
fprime = lambda x : np.tan(x) + x*(1/(np.cos(x)))**2

### Find initial guesses 
xrt_all = ()
for(guess in range):
  it, xroot, f_xroot = proot(guess, n, tol)
  xrt_all = xrt_all + (xroot,)
  
  

### Setup looped solution
range = np.arange(0,0.1,40)
tol = 0.000000001
n = 50


  
print "iterations: " + np.str(it)
#### prints the following:
#### "iterations: 6"
print "first root lower than: " + np.str(guess) + " is " + np.str(xroot)
#### prints the following:
#### "first root lower than: 1.4 is 1.19245882934"
print "at this root, the function evaluates as: " + np.str(f_xroot)
#### prints the following:
#### "at this root, the function evaluates as: -4.4408920985e-16"
