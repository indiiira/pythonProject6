
import math

import self as self

n = 10
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

def f(i):
  return f(float(i) * h)

def a(i):
  return p(float(i) * h)

def g(i):
  return q(float(i) * h)



  def __getitem__(self, i):
    return self.v[i]


def InitializeMatrix(alphas, betas):
   for i in range(n):
     for j in range(n):
        self.mtx[i][j] = 0.0
   ai = 0.0
   ci = a(0) + a(1) + h * h * g(0)
   bi = a(1)
   if (n > 2):
     alphas[0] = bi / ci
     betas[0] = 0.0
     mtx[0][0] = ai
     mtx[0][1] = ci
     mtx[0][2] = bi
   for i in range(1, n - 1):
     ai = self.alphas[i - 1] * a(i)
     ci = a(i) + a(i + 1) + h * h * g(i)
     bi = a(i + 1)
     mtx[i][i - 1] = ai
     mtx[i][i] = ci
     mtx[i][i + 1] = bi
     alphas[i] = bi / (ci - ai)
     betas[i] = (f(i) * h * h - self.betas[i - 1] * a(i)) / (ci - ai)
   ai = self.alphas[n - 2] * a(n - 1)
   ci = a(n - 1) + a(n) + h * h * g(n - 1)
   bi = 0.0
   alphas[n - 1] = bi / (ci - ai)
   betas[n - 1] = (f(n - 1) * h * h - self.betas[n - 2] * a(n - 1)) / (ci - ai)
   mtx[n - 1][n - 3] = ai
   mtx[n - 1][n - 2] = ci
   mtx[n - 1][n - 1] = bi


def getSolutions():
    self.y[n - 1] = 0.0
    self.y[0] = 0.0
    for i in range(n - 2, 0, -1):
      self.y[i] = self.alphas[i] * self.y[i + 1] + self.betas[i]
    return self.y

mtx = method1.getMatrix()
y1 = method1.getSolutions()
print("ih \t\t yi \t\t u(ih) \t\t |yi-u(ih)|")
for i in range(n):
 print(i * h, "\t\t", y1[i], "\t\t", u(i * h), "\t\t", abs(y1[i] - u(i * h)))
print("\n\n")