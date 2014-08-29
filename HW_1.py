# init

import numpy as np
import matplotlib.pyplot as plt
file = open("HW_1_output.txt","w")

###### Prob 1
file.write("Prob 1 Results:"+"\n")

### (a) new matrix with new row
file.write("(a):"+"\n")

A = np.array([[1,2,3,1],[4,1,3,2],[4,4,3,1],[3,2,3,3]])
z = np.array([[1,1,3,1]])
B = np.vstack((A[0:2,:],z,A[2:5,:]))

file.write(str(B)+"\n")

### (b) new matrix with new column
file.write("(b):"+"\n")
A = np.array([[1,2,3,1],[4,1,3,2],[4,4,3,1],[3,2,3,3]])
z = np.array([[1,1,3,1]])
C = np.hstack((A[:,0:2],z.T,A[:,2:5]))

file.write(str(C)+"\n"+"\n")

###### Prob 2
file.write("Prob 2 Results:"+"\n"+"\n")

import math as math

### (a) range
file.write("(a):"+"\n")
x = np.arange(-1*math.pi,math.pi,0.05*math.pi)

file.write(str(x)+"\n")

### (b) sin function
file.write("(b):"+"\n")
y = np.sin(x)
file.write(str(y)+"\n")

### (c) plot
file.write("(c):"+"Plot below"+"\n"+"\n")


#plt.plot(x,y)
#plt.show()

## THIS GENERATES A PLOT THAT I CAN ATTACH

###### Prob 3
file.write("Prob 3 Results:"+"\n")

a = np.ones(6)
b = np.ones(4)*-.05
M = np.diag(a)+np.diag(b,-2)

file.write(str(M)+"\n"+"\n")

###### Prob 4
file.write("Prob 4 Results:"+"\n")
rand = np.random.rand(5,5)

det = np.linalg.det(rand)

file.write(str(rand)+"\n")
file.write(str(det)+"\n"+"\n")

###### Prob 5
file.write("Prob 5 Results:"+"\n")

coef = np.array([[4,2,1],[3,1,1],[1,2,1]])
rhs = np.array([[2],[1],[1]])

soln = np.linalg.solve(coef, rhs)

file.write("x1:" + str(soln[0])[2:8] + " x2:" + str(soln[1])[2:8] + " x3:" + str(soln[2])[2:8]+"\n")

file.close()
