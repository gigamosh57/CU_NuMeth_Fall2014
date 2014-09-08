import numpy as np
import matplotlib.pyplot as plt
file = open("HW_2_output.txt","w")
## Problem 1
file.write("Prob 1 Results:"+"\n")
## if A is an [m x n] matrix, then B must be [d x n]  A/B is essentially solving a system of linear equations x * B = A
# Example calc
A = np.array([1,2,3])
B = np.array([[1,2,3],[2,4,8]])
x = np.divide(A,B)
file.write("np.dot(x,A): "+ str(np.dot(x,A))+"\n"+"\n")

## Problem 2
file.write("Prob 2 Results:"+"\n")
def mytri(n):
	a = np.diag(np.ones(n))
	b = np.diag(np.ones(n-1)*-0.5,-1)
	c = np.diag(np.ones(n-1)*-.25,1)
	A = a+b+c
	r = np.reshape(np.hstack((1,np.zeros(n-1))),(1,n)).T
	y = np.dot(np.linalg.inv(A),r)
	return A, y
A, y = mytri(5)
print A
file.write("A: " + str(A)+"\n")
print y
file.write("y: " + str(y)+"\n"+"\n")

## Problem 3
file.write("Prob 3 Results:"+"\n")
def mytri2(n,alpha):
	a = np.diag(np.ones(n))
	b = np.diag(np.ones(n-1)*-0.5,-1)
	c = np.diag(np.ones(n-1)*-.25,1)
	A = a+b+c
	r = np.reshape(np.hstack((alpha,np.zeros(n-1))),(1,n)).T
	y = np.dot(np.linalg.inv(A),r)
	return A, y
A, y = mytri2(5,2)
print A
file.write("A: " + str(A)+"\n")
print y
file.write("y: " + str(y)+"\n"+"\n")

## Problem 4
file.write("Prob 4 Results:"+"\n")
import math
##### NEED SCIPY FOR THIS TO WORK
##### INSTALL ANACONDA HERE
##### http://continuum.io/downloads
from scipy.optimize import fsolve
taneqn = lambda x : x * np.tan(x)-3
bounds =  np.mean([1.5,1])
ans = fsolve(taneqn,bounds)
file.write("ans: " + str(ans) + "\n"+"\n")

## Problem 5
import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt
# for fsolve
from scipy.optimize import fsolve
#file.write("Prob 5 Results:"+"\n")
#r = np.arange(0,10.5,0.5)
#tanroot = np.zeros(len(r))
taneqn2 = lambda x : x * np.tan(x) - r
#bd = np.mean([math.pi/2,0])
bd = np.mean([1.57,0])
#bd=0
tanroot=()
tanrootu=()
tanrootl=()
rplot=()
for r in np.arange(0,10.5,0.5):
    tanroot = tanroot + (fsolve(taneqn2,x0=bd,maxfev=1000),)
    tanrootu = tanrootu + (fsolve(taneqn2,x0=0,maxfev=1000),)
    tanrootl = tanrootl + (fsolve(taneqn2,x0=1.57,maxfev=1000),)
    rplot = rplot + (r,)

plt.plot(rplot,tanroot)
plt.plot(rplot,tanrootl)
plt.plot(rplot,tanrootu)
plt.show()
