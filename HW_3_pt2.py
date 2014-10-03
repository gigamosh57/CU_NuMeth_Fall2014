#### Page Weil
#### 10/1/14
#### CVEN 5537
#### HOMEWORK 3, PROBLEM 4

import numpy as np
# For pi
import math
# for erf
import scipy.special as sp
# for plotting
import matplotlib.pyplot as plt

################################################
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
################################################

# Define f
f_eb = lambda x : x*t - x**3

# Input this function into the EB method.  Every timestep, you need to solve:
#  x^n+1 - f^n+1*dt - x^n = 0
# converting this to a form that can be input into the Newton solver:
# f
f = lambda x : x - (x*t - x**3)*dt - xold
# fprime
fprime = lambda x : 1 - (t - 3*x**2)*dt

# Define timesteps
iter = 40
dt = 0.1
steps = np.arange(0,dt*iter,dt) #should be 40 steps

x0 = 1
xold = x0
xnews = ()
for t in steps:
  guess = xold
  xnew = mynewton(guess)
  # SOLVE fstep for xnew
  it, xnew, f_xroot = mynewton(guess)
  # Store xnew in list
  xnews = xnews + (xnew,)
  xold = xnew

# Using Wolfram Alpha to determine the analytical solution we find:
# x = np.exp(t**2/2)/np.sqrt(np.sqrt(np.pi*np.erfi(t)+1))

t_an = steps
x_an = np.exp(t_an**2/2)/np.sqrt(np.sqrt(np.pi)*sp.erfi(t_an)+1)

# Plot results

fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
axes.set_xlabel('t')
axes.set_ylabel('x(t), x0 = 1')

axes.plot(steps,x_an,color = "blue",label = "analytical")
axes.plot(steps,list(xnews),color = "red",label = "numerical")

axes.legend()
plt.show()



