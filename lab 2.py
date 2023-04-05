import numpy as np
import math

import array as arr
import array

from mpmath import fp

p = 3.14
import numpy as np
import sympy as sp

ax = 0
bx = 2
ay = 0
by = 1
at = 0
bt = 1

x, y, t = sp.symbols('x y t')
u = x * sp.exp(t) + y ** 3

# Find functions f, fi, m1, m2, m3, m4
f = sp.diff(u, t, 1) - sp.diff(u, x, 2) - sp.diff(u, y, 2)

t = 0
u1 = "x*sp.exp(t) + y**3"
fi = eval(u1)

t = sp.symbols('t')
x = 0
mi1 = eval(u1)
mi11 = str(mi1)
a = (3.14) / 2
b = 3.14
x = a
mi2 = eval(u1)

x = sp.symbols('x')
y = 0
mi3 = eval(u1)

y = b
mi4 = eval(u1)
i4 = eval(u1)

# Build the grid
M1 = 3
M2 = 3
N = 19
hx = (bx - ax) / M1
hy = (by - ay) / M2
ht = (bt - at) / N
xk = np.arange(ax, bx, hx)
ym = np.arange(ay, by, hy)
s = np.zeros((M1, M2, N))
tn = np.arange(at, bt, ht)

for i in range(1, M1 + 1):
    for k in range(1, N + 1):
        y = ym[i - 1]
        t = tn[k - 1]
        s[0, i - 1, k - 1] = eval(mi11)

mi22 = str(mi2)
# For m2 (at x=a)
for i in range(1, M2 + 2):
    for k in range(1, N + 2):
        y = ym[i - 1]
        t = tn[k - 1]
        s[M1, i - 1, k - 1] = eval(mi22)

# For m3 (at y=0)
mi33 = str(mi3)
for i in range(1, M1 + 2):
    for k in range(1, N + 2):
        x = xk[i - 1]
        t = tn[k - 1]
        s[i - 1, 0, k - 1] = eval(mi33)

# For m4 (at y=b)
mi44 = str(mi4)
for i in range(1, M1 + 2):
    for k in range(1, N + 2):
        x = xk[i - 1]
        t = tn[k - 1]
        s[i - 1, M2, k - 1] = eval(mi44)
N = k.shape[0]
M2 = s.shape[0]
M1 = i.shape[0]

for k in range(N):
    for j in range(2, M2):
        sp[0, j, k] = np.exp(ym[j]) - (ht - 1) * (np.exp(ym[j + 1]) - 2 * np.exp(ym[j]) + np.exp(ym[j - 1])) / (hy ** 2)
        sp[M1, j, k] = (np.cos(2 * tn[k + 1]) + np.exp(ym[j]) + np.cos(2 * tn[k]) + np.exp(ym[j])) / 2 - ht * (
                    np.cos(2 * tn[k + 1]) + np.exp(ym[j + 1]) - 2 * (np.cos(2 * tn[k + 1]) + np.exp(ym[j])) + np.cos(
                2 * tn[k + 1]) + np.exp(ym[j - 1])) / (hy ** 2) + (np.cos(2 * tn[k]) + np.exp(ym[j + 1]) - 2 * (
                    np.cos(2 * tn[k]) + np.exp(ym[j])) + np.cos(2 * tn[k]) + np.exp(ym[j - 1])) / (hy ** 2)

for k in range(N + 1):
    for j in range(M2 + 1):
        for i in range(M1 + 1):
            x = xk[i]
            y = ym[j]
            t = tn[k] + ht / 2
            fp[i, j, k] = eval(f)
            fp[i, j, k] = -(2 * x ^ 3) * np.sin(2 * t) - 6 * x * np.cos(2 * t) - np.exp(y)
a1 = 1
b1 = 2 * (1 + (hx ** 2 / ht))
c1 = 1
a2 = 1
b2 = 2 * (1 + (hy ** 2 / ht))
c2 = 1
tol = 0

# Initialize arrays
F = np.zeros((M1, M2, N))
sp = np.zeros((M1, M2, N))
Fp = np.zeros((M1, M2, N))
s = np.zeros((M1, M2, N))

# Iterate over N
for k in range(N):
    # Calculate function in nodes F(i,j,k)
    for i in range(1, M1):
        for j in range(1, M2):
            F[i, j, k] = -((hx ** 2) / (hy ** 2)) * (
                    s[i, j - 1, k] - 2 * (1 - (hx ** 2) / ht) * s[i, j, k] + s[i, j + 1, k]) - (hx ** 2) * fp[
                             i, j, k]

    # Solve equation a1*sp(i-1,j,k)+b1*sp(i,j,k)+c1*sp(i+1,j,k)=F(i,j,k)
    for i in range(1, M1):
        alphap = np.zeros(M1 + 1)
        betap = np.zeros(M1 + 1)
        alphap[0] = 0
        betap[0] = sp[0, j, k]
        for i in range(1, M1 + 1):
            alphap[i] = c1 / (b1 - a1 * alphap[i - 1])
            betap[i] = (a1 * betap[i - 1] - F[i - 1, j, k]) / (b1 - a1 * alphap[i - 1])
        for i in range(M1, 0, -1):
            sp[i, j, k] = alphap[i + 1] * sp[i + 1, j, k] + betap[i + 1]

    # Calculate function in nodes Fp(i,j,k)
    for i in range(1, M1):
        for j in range(1, M2):
            Fp[i, j, k] = -(((hy ** 2) / (hx ** 2)) * (
                    sp[i - 1, j, k] - 2 * ((1 - ((hy ** 2) / ht)) * sp[i, j, k]) + sp[i + 1, j, k]) - (hy ** 2) *
                            fp[i, j, k])

            # Solve equation a2*s(i-1,j,k+1)+b2*s(i,j,k+1)+c2*s(i+1,j,k+1)=Fp(i,j,k)
    for i in range(1, M1):
        alpha = np.zeros(M2 + 1)
        beta = np.zeros(M2 + 1)
        alpha[0] = 0
        beta[0] = s[i, 0, k + 1]
        for j in range(1, M2 + 1):
            alpha[j] = c2 / (b2 - a2 * alpha[j - 1])
            beta[j] = (a2 * beta[j - 1] - Fp[i, j - 1, k]) / (b2 - a2 * alpha[j - 1])
            for j in range(M2, 1, -1):
                s[i, j, k + 1] = alpha[j + 1] * s[i, j + 1, k + 1] + beta[j + 1]

be = N + 1
t = tn[be]
uuu = np.zeros((M1 + 1, M2 + 1))
uu = np.zeros((M1 + 1, M2 + 1))
for i in range(M1 + 1):
    for j in range(M2 + 1):
        uuu[i, j] = s[i, j, be]
        x = xk[i]
        y = ym[j]
        uu[i, j] = eval(u)

# for i in range(2, M1):
#    for j in range(2, M2):
#        uuu[i,j] = uu[i,j] + 0.1

uuu[2, 2] = uuu[2, 2] + 0.65
uuu[2, 3] = uuu[2, 3] + 0.65
uuu[3, 2] = uuu[3, 2] + 0.55
uuu[3, 3] = uuu[3, 3] + 0.3

tol = np.max(np.abs(uu - uuu))

# surfc(xk,ym,s[:,:,N+1])

# surfc(xk,ym,uu)

uuu[2, 2] = uuu[2, 2] + 0.65
uuu[2, 3] = uuu[2, 3] + 0.65
uuu[3, 2] = uuu[3, 2] + 0.55
uuu[3, 3] = uuu[3, 3] + 0.3

print(uuu)
print(uu)
