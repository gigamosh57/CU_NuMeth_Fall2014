


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
n = 50
it, xroot, f_xroot = mynewton(guess, n, tol)
print "iterations: " + np.str(it)
#### prints the following:
#### "iterations: 6"
print "first root lower than: " + np.str(guess) + " is " + np.str(xroot)
#### prints the following:
#### "first root lower than: 2 is 0.48539431657"
print "at this root, the function evaluates as: " + np.str(f_xroot)
#### prints the following:
#### "at this root, the function evaluates as: -8.32667268469e-17"


# Define second Function
f = lambda x : x*np.tan(x)-3
# Fprime
fprime = lambda x : np.tan(x) + x*(1/(np.cos(x)))**2

### Solve with initial guess
guess = 1.4
tol = 0.000000001
n = 50
it, xroot, f_xroot = mynewton(guess, n, tol)
print "iterations: " + np.str(it)
#### prints the following:
#### "iterations: 6"
print "first root lower than: " + np.str(guess) + " is " + np.str(xroot)
#### prints the following:
#### "first root lower than: 1.4 is 1.19245882934"
print "at this root, the function evaluates as: " + np.str(f_xroot)
#### prints the following:
#### "at this root, the function evaluates as: -4.4408920985e-16"

