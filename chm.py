import math
import numpy as np
from tabulate import tabulate


def Px(x):
    return 1+x**2

def Qx(x):
    return x + 1

def Ux(x):
    return x**4*((1-x)**2)

def Fx(x):
    return -x**7-43*x**6-59*x**5-49*x**4-40*x**3-12*x**2


n = 10
h = 1.0 / n
eps = h ** 6

x = np.zeros(n+1)
ai = np.zeros(n+1)
qi = np.zeros(n+1)
fi = np.zeros(n+1)
alpha = np.zeros(n+1)
beta = np.zeros(n+1)
a_big = np.zeros(n+1)
b_big = np.zeros(n+1)
c_big = np.zeros(n+1)
A = np.zeros((3, n+1))
y = np.zeros(n+1)
y_jacobi = np.zeros(n+1)
y_relax = np.zeros(n+1)
alpha[0] = 0
beta[0] = 0
y[0] = 0
y[n] = 0
k = 0
best_k = 99999
best_omega = 0
for i in range(n+1):
    ai[i] = Px(i * h)
    qi[i] = Qx(i * h)
    fi[i] = Fx(i * h) * h ** 2
    x[i] = i * h

for i in range(1, n):
    a_big[i] = -ai[i]
    b_big[i] = -(ai[i] + ai[i + 1] + h ** 2 * qi[i])
    c_big[i] = -ai[i + 1]

for i in range(n):
    A[0, i] = a_big[i]
    A[1, i] = b_big[i]
    A[2, i] = c_big[i]

for i in range(1, n):
    alpha[i+1] = c_big[i] / (b_big[i] - a_big[i] * alpha[i])
    beta[i+1] = (a_big[i] * beta[i] - fi[i]) / (b_big[i] - a_big[i] * alpha[i])

for i in range(n-1, 0, -1):
    y[i] = alpha[i+1] * y[i+1] + beta[i+1]

def Jacobi(a, f, g):
    global k
    y = np.zeros(n+1)
    y_pre = np.zeros(n+1)
    r = 0

    for i in range(n):
        y[i] = 0

    while r > eps or k == 0:
        for i in range(n):
            y_pre[i] = y[i]

        r = -1

        for i in range(1, n):
            y[i] = (f[i]*h**2+ a[i+1] * y_pre[i+1] + a[i] * y_pre[i-1]) / (a[i] + a[i+1] + g[i] * h ** 2)

        k += 1

        for i in range(1, n):
            if abs(-a[i+1] * y[i+1] + (a[i+1] + a[i] + g[i] * h ** 2) * y[i] - a[i] * y[i-1] - f[i]*h**2) > eps:
                r = abs(-a[i+1] * y[i+1] + (a[i+1] + a[i] + g[i] * h ** 2) * y[i] - a[i] * y[i-1] - f[i]*h**2)
    return y
y_relax = Jacobi(ai, fi, qi)


print(f"Вычисления для омега={best_omega}\n")
for i in range(n+1):
        print(tabulate([[i * h, y[i],y_relax[i], abs(y[i] - y_relax[i])]], headers=['i*h', 'y[i]', 'U(i*h)', '|y[i]-U(i*h)|']))



