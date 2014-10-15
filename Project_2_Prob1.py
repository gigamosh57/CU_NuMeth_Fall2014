#### Page Weil
#### 10/8/14
#### CVEN 5537
#### PROJECT 2 , PROBLEM 1

import numpy as np
# For pi
import math
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

################################################
######## RK4 ODE solver begins here
###### I modified the codes provided so they 
###### would work in Python

#general-purpose function for solving a system of ODEs using RK4 method
#you need to provide:
# f - function file name that contains the derivative functions
#(note that f should evaluate the derivative vector as a column vector)
#x0 - initial condition (column vector); tstart,tfinal - start and end times, dt - time step
#the output is a time vector and a matrix xsol - each row of xsol contains
#the transpose of the x vector at the corresponding time

def myrk4(f,x0,tstart,tfinal,dt):
  #initialize time vector
  tvec=[tstart]
  #initialize x, initialize solution matrix xsol
  x=x0
  xsol=[] + [x,]
  time = np.arange(tstart,tfinal,dt)
  for t in time:
    #stage 1
    k1=f(t,x)*dt;
    #stage 2
    k2=f(t+dt/2,x+k1/2)*dt;
    #stage 3
    k3=f(t+dt/2,x+k2/2)*dt;
    #stage 4
    k4=f(t+dt,x+k3)*dt;
    #calculate x at next time-step (using same variable name x, I don't need xold, xnew etc.)
    x=x+(k1+2*k2+2*k3+k4)/6;
    #update time vector and xsol matrix
    tvec=tvec + [t+dt]
    xsol=xsol + [x]
  return tvec,xsol
######## RK4 ODE solver ends here
################################################


#### Begin Problem 1

dts = [0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5]#np.arange(0.01,1,0.01)
gte1 = []
gte2 = []

for a in dts:
  f = lambda t,x: -x  
  tvec, xsol = myrk4(f,1,0,1,a)
  # This is the analytical solution for f = -x
  fana = lambda t: np.exp(-1*np.array(t))
  gte1 = gte1 + [xsol[len(xsol)-1] - fana(tvec[len(xsol)-1])]
  
  f = lambda t,x: -x**2  
  tvec, xsol = myrk4(f,1,0,1,a)
  # This is the analytical solution for f = -x^2
  fana = lambda t: 1/(1+np.array(t))
  #gte2 = gte2 + [sum(xsol - fana(tvec))]
  gte2 = gte2 + [xsol[len(xsol)-1] - fana(tvec[len(xsol)-1])]

fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
px = np.log(gte1)
py = np.log(dts)

axes.plot(px,py,'ro',label="f = -x")
fit = np.polyfit(px,py,1)
axes.plot(px,fit[0]*px+fit[1],label="slope = " + str(1/fit[0]))

px = np.log(gte2)
py = np.log(dts)
axes.plot(px,py,'bo',label="f = -x^2")
fit = np.polyfit(px,py,1)
axes.plot(px,fit[0]*px+fit[1],label="slope = " + str(1/fit[0]))

axes.set_xlabel('ln(GTE)')
axes.set_ylabel('ln(delta t)')
axes.legend()
plt.show()

#Based on the slope in the picture shown, we have found that the error of RK4 is o(delta t ^5)

