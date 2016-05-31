######################################################################################
#    Program to solve the Eikonal equation in 2D using the Rouy and Tourin Scheme
#         u_t + |\Grad u_x| = 1 in ]-1,1[ x
#           u(x,0) = 0
#           u = 0 on boundary
#
#
#
######################################################################################

import sys
sys.modules[__name__].__dict__.clear()
import numpy as np
import matplotlib.pylab as plt
import time

a = -1.
b = 1.
c = -1.
d = 1.

N = 100

h = (b-a)/N

x = np.zeros(N-1)
y = np.zeros(N-1)

for j in range(N-1):
    x[j] = a + (j+1)*h
    y[j] = c + (j+1)*h

#CFL condition
delt = 0.1*h

u = np.zeros((N-1,N-1))
unew = np.zeros((N-1,N-1))

iterations = 0
error = 100
tol = 1e-5

fig = plt.figure()
t0 = time.time()
while error > tol:
    for i in range(1,N-2):
        for j in range(1,N-2):
            Dxp = (u[i+1,j]-u[i,j])/h
            Dxm = (u[i-1,j]-u[i,j])/h
            Dyp = (u[i,j+1]-u[i,j])/h
            Dym = (u[i,j-1]-u[i,j])/h

            Dx = np.min([0,Dxp,Dxm])
            Dy = np.min([0,Dym,Dyp])

            H = np.sqrt(Dx**2+Dy**2)-1

            unew[i,j] = u[i,j] - delt*H

    error = np.amax(abs(u-unew))
    np.copyto(u,unew)

    iterations = iterations + 1

t1 = time.time()

print 'Total Time : ', t1-t0
plt.contour(x,y,u)
plt.draw()
