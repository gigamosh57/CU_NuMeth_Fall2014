#### NUMERICAL METHODS PROJECT 1, PROBLEM 3
import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt

dxdt = lambda x : -x
x0 = 1
dt = 0.9
t = np.arange(0,20,dt)
xnow = x0
for a in t:
	xnext = xnow + dxdt(xnow)*dt
	xnow = xnext
	print xnext
