############################################################
#       Code to solve the advection equation in 1D
#           u_t + u_x = 0 in ]-10,10[ x (0,T]
#               u(x,0) = u0(x) = sin(pi*x)
#           u(-10,t) = u(10,t) = 0
#
#   Exact solution : u(x,t) = u0(x-t) = sin(pi*(x-t))
###########################################################

import sys
sys.modules[__name__].__dict__.clear()

import numpy as np
import matplotlib.pylab as plt


a = -1.
b = 4.

N = 100 # N subintervals

h = (b-a)/N

#CFL Condition
delt = 0.8*h

tf = 1
t0 = 0

M = (tf-t0)/delt

x = np.zeros(N-1)
for j in range(N-1):
    x[j] = a + (j+1)*h

u = np.zeros(N-1)
unew = np.zeros(N-1)
uexact = np.zeros(N-1)

lamb = delt/h;

for i in range(N-1):
    if x[i] < 0:
        u[i] = 1
    elif x[i] >= 0 and x[i] <= 1:
        u[i] = 2*x[i]**3 - 3*x[i]**2 + 1
    else:
        u[i] = 0

plt.ion()
for k in range(int(M)):
    unew[0] = u[0] - lamb*(u[0]-1)
    for j in range(1,N-1):
        unew[j] = u[j] - lamb*(u[j]-u[j-1])

    for i in range(N-1):
        u[i] = unew[i]

    plt.plot(x,unew)
    plt.draw()
    plt.pause(0.2)

for i in range(N-1):
    if x[i]-tf < 0:
        uexact[i] = 1
    elif x[i]-tf >= 0 and x[i]-tf <= 1:
        uexact[i] = 2*(x[i]-tf)**3 - 3*(x[i]-tf)**2 + 1
    else:
        uexact[i] = 0

plt.plot(x,unew,'g',x,uexact,'r')
plt.draw()
