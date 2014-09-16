#### NUMERICAL METHODS PROJECT 1, PROBLEM 5
import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt

T_0 = 0
T_d = 1
L = 1
alpha = 1
nmax = 50
theta0 = 1

theta = lambda x, n, t : theta0*(1-4*(np.sin(((2*n+1)*np.pi*x)/(2*L))*np.exp(-(((2*n+1)**2)*(np.pi**2)*alpha*t)/(4*L**2)))/((2*n+1)*math.pi))




dt = 0.01
t = np.arange(0,20,dt)
n = np.arange(0,nmax+1)


