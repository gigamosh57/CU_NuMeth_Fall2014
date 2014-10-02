#### Page Weil
#### 10/1/14
#### CVEN 5537
#### HOMEWORK 3, PROBLEM 4

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt

######## Newton  solver begins here
def mynewton(guess,n = 1000,tol = 0.01): 
  it = 0
  i = 0
  xnow = guess
  while (i < n):
      i = i + 1
      it = i
      deltax = -1*f(xnow)/fprime(xnow)
      xnow = xnow + deltax
      if (np.absolute(f(xnow)) < tol):
            i = n+1
  return it, xnow, f(xnow)
######## Newton solver ends here

# Define First Function
f = lambda x : x*t - x**3
# Fprime
fstep = lambda x : x*t - x**3 - xold

### Solve with initial guess
guess = 1.4
tol = 0.000000001
n = 50
it, xroot, f_xroot = mynewton(guess, n, tol)

iter = 40
tstep = 0.1
steps = np.arange(tstep,iter) #should be 40 steps

x0 = 1
xold = x0
for a in steps:
  # SOLVE fstep for xnew
  # Store xnew in list
  # xold = xnew
  
# plot results




