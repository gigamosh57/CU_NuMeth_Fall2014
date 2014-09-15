#### NUMERICAL METHODS PROJECT 1, PROBLEM 3
import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt

######## Bisection solver begins here
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
  
######## Bisection solver ends here



# Define Function
f = lambda x : x*np.tan(x)-3
# Fprime
fprime = lambda x : np.tan(x) + x*(1/(np.cos(x)))**2

### Find initial guesses 
inc = 0.1
rng = np.arange(0,40,inc)
a = zip(rng,rng+inc)
n = 1000
tol = 0.0000001

### Loop through range and find initial guesses:
guessinit = ()
for g in a:
#  print g
  guess = list(g)
  t = [f(guess[0]),f(guess[1])] 
  if (all( item > 0 for item in t) or all( item <= 0 for item in t)) == False:
    guessinit = guessinit + (guess,)

### Loop through guesses and find roots
xrt_all = ()
fxrt = ()
for g in guessinit:
  guess = list(g)
  it, xroot, f_xroot = proot(guess, n, tol)
  xrt_all = xrt_all + (xroot,)
  fxrt = fxrt + (f_xroot,)
  
  
### Loop through results and find actual roots (not +/- infinity errors)


