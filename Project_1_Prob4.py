#### NUMERICAL METHODS PROJECT 1, PROBLEM 3
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
xplot = np.arange(0,len(v),1)/float(len(v))
plt.plot(xplot,v,color = "b")

v = x_05	
xplot = np.arange(0,len(v),1)/float(len(v))
plt.plot(xplot,v,color = "g")

v = x_1
xplot = np.arange(0,len(v),1)/float(len(v))
plt.plot(xplot,v,color = "r")

v = x_2
xplot = np.arange(0,len(v),1)/float(len(v))
plt.plot(xplot,v,color = "c")

v = x_25
xplot = np.arange(0,len(v),1)/float(len(v))
plt.plot(xplot,v,color = "m")

plt.show()

### This plot show s



