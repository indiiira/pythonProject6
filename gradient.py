

import math
import numpy as np
from tabulate import tabulate

n = 10
h = 1.0 / n
eps = h ** 3

x = np.zeros(n+1)
ai = np.zeros(n+1)
qi = np.zeros(n+1)
fi = np.zeros(n+1)
Alpha = np.zeros(n+1)
Beta = np.zeros(n+1)
a_big = np.zeros(n+1)
b_big = np.zeros(n+1)
c_big = np.zeros(n+1)
A = np.zeros((3, n+1))
y = np.zeros(n+1)
y_jacobi = np.zeros(n+1)
y_relax = np.zeros(n+1)

alpha = 4.0
beta = 2.0
gamma = 2.0
h = 1.0 / float(n)
EPS = h * h * h

def p(x):
  return 1 + math.pow(x, gamma)

def p_(x):
  return gamma * (1 / math.pow(x, gamma - 1))

def q(x):
  return x + 1

def u(x):
  return math.pow(x, alpha) * math.pow(1 - x, beta)

def u_(x):
  return alpha * beta * (math.pow(x, alpha - 1)) * (math.pow(1 - x, beta - 1)) - beta * math.pow(x, alpha) * (math.pow(1 - x, beta - 1))

def u__(x):
  return alpha * (alpha - 1) * (math.pow(x, alpha - 2)) * math.pow(1 - x, beta) - 2 * alpha * beta * (math.pow(x, alpha - 1)) * (math.pow(1 - x, beta - 1)) + beta * (beta - 1) * math.pow(x, alpha) * (math.pow(1 - x, beta - 2))

def f(x):
  return -(p_(x) * u_(x) + p(x) * u__(x)) + q(x) * u(x)
k = 0
best_k = 99999
best_omega = 0

def Px(x):
    return 1 + x ** 2

def Qx(x):
    return x + 1

def Ux(x):
    return x**4*((1-x)**2)

def Fx(x):
    return -x**7-43*x**6-59*x**5-49*x**4-40*x**3-12*x**2

Alpha[0] = 0
Beta[0] = 0
y[0] = 0
y[n] = 0

for i in range(1,n+1):
    ai[i] = p(i * h)
    qi[i] = q(i * h)
    fi[i] = f(i * h) * h ** 2
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
    Alpha[i+1] = c_big[i] / (b_big[i] - a_big[i] * Alpha[i])
    Beta[i+1] = (a_big[i] * Beta[i] - fi[i]) / (b_big[i] - a_big[i] * Alpha[i])

for i in range(n-1, 0, -1):
    y[i] = Alpha[i+1] * y[i+1] + Beta[i+1]

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
            y[i] = (f[i] + a[i+1] * y_pre[i+1] + a[i] * y_pre[i-1]) / (a[i] + a[i+1] + g[i] * h ** 2)

            k += 1

        for i in range(1, n):
            if abs(-a[i+1] * y[i+1] + (a[i+1] + a[i] + g[i] * h ** 2) * y[i] - a[i] * y[i-1] - f[i]) > r:
                r = abs(-a[i+1] * y[i+1] + (a[i+1] + a[i] + g[i] * h ** 2) * y[i] - a[i] * y[i-1] - f[i])

    return y

y_jacobi = Jacobi(ai, fi, qi)

with open("yakobi.txt", "w") as f:
    for i in range(n+1):
        f.write(f"{i*h}\t{y[i]}\t{y_jacobi[i]}\t{abs(y[i]-y_jacobi[i]):.10f}\n")
    f.write(f"Количество итераций: {k}")

def Relax(a, f, g):
    global k
    global best_k
    global best_omega
    y = np.zeros(n+1)
    y_pre = np.zeros(n+1)
    r = 0

    for omega in np.arange(0.05, 1, 0.05):
        for i in range(n+1):
            y[i] = 0

        k = 0

        while r > eps or k == 0:
            for i in range(n+1):
                y_pre[i] = y[i]

            r = -1

            for i in range(1, n):
                y[i] = ((f[i] + a[i+1] * y_pre[i+1] + a[i] * y_pre[i-1]) / (a[i+1] + a[i] + g[i] * h ** 2)) * omega + (1 - omega) * y_pre[i]

            k += 1

            for i in range(1, n):
                if abs(-a[i+1] * y[i+1] + (a[i+1] + a[i] + g[i] * h ** 2) * y[i] - a[i] * y[i-1] - f[i]) > r:
                    r = abs(-a[i+1] * y[i+1] + (a[i+1] + a[i] + g[i] * h ** 2) * y[i] - a[i] * y[i-1] - f[i])

        if k < best_k:
            best_k = k
            best_omega = omega
            y_relax = np.copy(y)

    return y_relax

y_relax = Relax(ai, fi, qi)


with open("relax.txt", "w") as f:
    f.write(f"Вычисления для омега={best_omega}\n")
    for i in range(n+1):
        f.write(f"{i*h}\t{y[i]}\t{y_relax[i]}\t{abs(y[i]-y_relax[i]):.10f}\n")
    f.write(f"Количество итераций: {k}")
for i in range(len(y)):
       print(tabulate([[i * h, y[i], Ux(i * h), abs(y[i] - Ux(i * h))*h*h]],   headers=['i*h', 'прогонка', 'u(x)', '|y[i]-U(i*h)|']))



def S(a, f, g):

    global k
    y = np.zeros(n + 1)
    y_pre = np.zeros(n + 1)
    r = np.zeros(n + 1)
    Ar=0
    kr=0
    for i in range(1, n):
      r[i] = -a[i] * y[i - 1] + a[i] * y[i] + a[i + 1] * y[i] + (h ** 2) * g[i] * y[i] - a[i + 1] * y[i + 1] - f[i] * h ** 2

      if abs(r[i])<eps:
          y[i] = y[i - 1] - Ar * r[i]
          k = k + 1
      else:
          for i in range(1, n):
              Ar = -a[i] * r[i - 1] + (a[i] + a[i + 1] + (h ** 2) * g[i]) * r[i] - a[i + 1] * r[i + 1]
              y[i] = y[i - 1] - Ar * r[i]
              k = k + 1



    return y

y_spusk=S(ai, fi, qi)

for i in range(n+1):
        print(tabulate([[i * h, y[i], y_spusk[i], abs(y[i] - y_spusk[i])]], headers=['i*h', 'y[i]', 'spusk', '|y[i]-U(i*h)|']))
print(f"Количесвто итераций {k}")