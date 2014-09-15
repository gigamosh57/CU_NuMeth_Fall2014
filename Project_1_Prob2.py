


## This code calculates the roots of an equation from prob 5, hw 2
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
      print xnow, f(xnow)
  return it, xnow, f(xnow)
######## Newton solver ends here



# Define First Function
f = lambda x : x**3 - np.cos(3*x)
# Fprime
fprime = lambda x : 3*x**2 + 3*np.sin(3*x)

### Solve with initial guess
guess = 2
tol = 0.000000001
n = 1000
it, xroot, f_xroot = mynewton(guess, n, tol)
print "iterations: " + np.str(it)
print "first root lower than: " + np.str(guess) + " is " + np.str(xroot)
print "at this root, the function evaluates as: " + np.str(f_xroot)


# Define second Function
f = lambda x : x*np.tan(x)-3
# Fprime
fprime = lambda x : np.tan(x) + x*np.arccos(x)

### Solve with initial guess
guess = 2
tol = 0.000000001
n = 1000
it, xroot, f_xroot = mynewton(guess, n, tol)
print "iterations: " + np.str(it)
print "first root lower than: " + np.str(guess) + " is " + np.str(xroot)
print "at this root, the function evaluates as: " + np.str(f_xroot)

