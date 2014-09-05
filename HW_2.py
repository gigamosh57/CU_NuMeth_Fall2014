import numpy as np

import matplotlib.pyplot as plt

#file = open("HW_1_output.txt","w")

## Problem 1
## if A is an [m x n] matrix, then B must be [d x n]  A/B is essentially solving a system of linear equations x * B = A
# Example calc
A = np.array([1,2,3])
B = np.array([[1,2,3],[2,4,8]])
x = np.divide(A,B)
print np.dot(x,A)

## Problem 2
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
print y
