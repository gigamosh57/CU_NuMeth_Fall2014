#### NUMERICAL METHODS PROJECT 1, PROBLEM 1
import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt

# Define Function
f = lambda x : x**3 - np.cos(3*x)
# Fprime
fprime = lambda x : 3*x**2 + 3*np.sin(3*x)

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

### TEST SOLUTION FOR SINGLE VALUE
guess = 2
tol = 0.000000001
n = 50
it, xroot, f_xroot = mynewton(guess, n, tol)

print "iterations: " + np.str(it)
print "first root lower than: " + np.str(guess) + " is " + np.str(xroot)
#### this prints
#### "first root lower than: 2 is 0.48539431657"
print "at this root, the function evaluates as: " + np.str(f_xroot)
#### this prints
#### "at this root, the function evaluates as: -8.32667268469e-17"

