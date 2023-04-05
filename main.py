import math

import numpy as np
from pip._internal.utils.misc import tabulate
from tabulate import tabulate
E = 2.7182818284959


def qx(x):
    return x+1


def px(x):
    return 1+x**2


def fx(x):
    return -x**7-43*x**6-59*x**5-49*x**4-40*x**3-12*x**2


def ux(x):
    return x**4*((1-x)**2)


def main():

    n = 10


    ai = np.zeros(n+1)
    qi = np.zeros(n+1)
    fi = np.zeros(n+1)
    alpha = np.zeros(n+1)
    beta = np.zeros(n+1)
    A = np.zeros(n + 1)
    B = np.zeros(n + 1)
    C = np.zeros(n + 1)
    y = np.zeros(n+1)

    h = 1.0 / n

    alpha[0] = 0
    beta[0] = 0
    y[0] = 0
    y[n] = 0


    # вычисление a[i],q[i],f[i]
    for i in range(n + 1):
        ai[i] = px(i * h)
        qi[i] = qx(i * h)
        fi[i] = fx(i * h) * h**2



    for i in range(1, n):
       A[i] = -ai[i]
       B[i] = -(ai[i] + ai[i + 1] + h * h * qi[i])
       C[i] = -ai[i + 1]

    # вычисление alpha, beta
    for i in range(1, n):
        alpha[i + 1] = -C[i] / (B[i] -A[i] * alpha[i])
        beta[i + 1] = (A[i] * beta[i]-fi[i]) / (B[i] - A[i] * alpha[i])

    # вычисление y[i]
    for i in range(n - 1, 0, -1):
        y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1]

    max = -100

    for i in range(n):
            if (math.fabs(y[i] -ux(i * h)) > max):
                max = math.fabs(y[i] - ux(i * h))



    for i in range(n + 1):

         print(tabulate([[i * h, y[i], ux(i * h),abs(y[i] - ux(i * h)*h) ]], headers=['i*h',	'y[i]',	'U(i*h)',	'|y[i]-U(i*h)|']))
    print(max)





if __name__ == "__main__":
    main()