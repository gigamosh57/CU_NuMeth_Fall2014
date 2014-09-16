#### NUMERICAL METHODS PROJECT 1, PROBLEM 4
import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt

dxdt = lambda x : -x
x0 = 1

dt = 0.1
t = np.arange(0,20,dt)
xnow = x0
x_01 = ()
for a in t:
	xnext = xnow + dxdt(xnow)*dt
	xnow = xnext
	x_01 = x_01 + (xnow,)

dt = 0.5
t = np.arange(0,20,dt)
xnow = x0
x_05 = ()
for a in t:
	xnext = xnow + dxdt(xnow)*dt
	xnow = xnext
	x_05 = x_05 + (xnow,)

dt = 1
t = np.arange(0,20,dt)
xnow = x0
x_1 = ()
for a in t:
	xnext = xnow + dxdt(xnow)*dt
	xnow = xnext
	x_1 = x_1 + (xnow,)

dt = 2
t = np.arange(0,20,dt)
xnow = x0
x_2 = ()
for a in t:
	xnext = xnow + dxdt(xnow)*dt
	xnow = xnext
	x_2 = x_2 + (xnow,)

dt = 2.5
t = np.arange(0,20,dt)
xnow = x0
x_25 = ()
for a in t:
	xnext = xnow + dxdt(xnow)*dt
	xnow = xnext
	x_25 = x_25 + (xnow,)

v = x_01	
xplot = np.arange(0,len(v),1)*20/float(len(v))
plt.plot(xplot,v,color = "b")

v = x_05	
xplot = np.arange(0,len(v),1)*20/float(len(v))
plt.plot(xplot,v,color = "g")

v = x_1
xplot = np.arange(0,len(v),1)*20/float(len(v))
plt.plot(xplot,v,color = "r")

v = x_2
xplot = np.arange(0,len(v),1)*20/float(len(v))
plt.plot(xplot,v,color = "c")

v = x_25
xplot = np.arange(0,len(v),1)*20/float(len(v))
plt.plot(xplot,v,color = "m")

plt.show()

### This plot shows the following behaviors:
# dt = 0.1, Solution converges to zero over time
# dt = 0.5, Solution converges to zero over time but converges to zero more quickly 
#  (this convergence is probably due to the larger steps)
# dt = 1, Result is always zero.  Since the initial condition is equal to 1, the 
#   function gets zero as its first answer and cannot diverge or converge further
# dt = 2, Result is an oscillation between -1 and 1.  This is because the convergence always 
#   attempts to move closer to zero but, because dt = x0 * 2 the function moves around zero 
#   back and forth to infinity
# dt = 2.5, Result is an oscillation positive and negative numbers that diverges over time.  
#   Every time the function moves closer to zero it moves farther and it moves past zero




