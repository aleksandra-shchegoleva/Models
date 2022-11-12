import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import *

mpl.rcParams.update({'font.size': 18})

N = 600
h = 0.1
M = np.arange(0, N, h)
x = [[1, 1, 1]]
h = 0.1
r, K, lambda1, lambda2, lambda3, mu = 0.1, 10, 0.04, 0.5, 0.35, 0.1
print(f"{mu*(lambda1 + lambda3)}, {sqrt(-mu*(lambda1 - lambda3)*(4*lambda1*lambda3 - lambda1*mu + lambda3*mu))}")
print(f'lambda1 > lambda3 : {lambda1}, {lambda3} или lambda1 > (-lambda3*mu)/(4*lambda3-mu) : {(-lambda3*mu)/(4*lambda3-mu)}\n lambda1 < lambda3 - mu : {lambda3 - mu}')
if (lambda1 > lambda3 or lambda1 > (-lambda3*mu)/(4*lambda3-mu)) and x[0][0] == x[0][1] and lambda1 < lambda3 - mu:
    x1c = -(mu*(lambda1 + lambda3) + sqrt(-mu*(lambda1 - lambda3)*(4*lambda1*lambda3-lambda1*mu+lambda3*mu)))/(2*(lambda1-lambda3+mu))
    x2c = x1c
#     x3c = r*(K*lambda1+K*x1c-x1c*(lambda1 + x1c))/(K*x1c)
    x3c = (K*lambda1*lambda2*r + K*lambda1*r*x1c - K*lambda1*x1c**2 + K*lambda2*r*x1c + K*r*x1c**2 - K*x1c**3 - lambda1*lambda2*r*x1c - lambda1*r*x1c**2 - lambda2*r*x1c**2 - r*x1c**3)/(K*x1c*(lambda2 + x1c))
    stpoint_1 = (-mu*(lambda1 + lambda3)/2 - sqrt(-mu*(lambda1 - lambda3)*(4*lambda1*lambda3 - lambda1*mu + lambda3*mu))/2)/(lambda1 - lambda3 + mu)
    stpoint_2 = (-mu*(lambda1 + lambda3)/2 + sqrt(-mu*(lambda1 - lambda3)*(4*lambda1*lambda3 - lambda1*mu + lambda3*mu))/2)/(lambda1 - lambda3 + mu)
elif lambda1 > lambda3 - mu:
    x1c = 0
    x2c = K
    x3c = 0
#     print(stpoint_1, stpoint_2)

text = f"$x_{1}(0)$ = {x[0][0]}\n" \
       f"$x_{2}(0)$ = {x[0][1]}\n" \
       f"$x_{3}(0)$ = {x[0][2]}\n" \
       f"$r$ = {r}\n" \
       f"$K$ = {K}\n" \
       f"$\mu$ = {mu},\n" \
       f"$\lambda_{1}$ = {lambda1}\n" \
       f"$\lambda_{3}$ = {lambda3}\n" \
       f"$h$" + f" = {h}"

for i in range(len(M) - 1):
    x.append([0, 0, 0])
    bx1, bx2 = r*(1 - x[i][0]/K), r*(1 - x[i][1]/K)
    fx1, fx2 = x[i][0]/(lambda1 + x[i][0]), x[i][1]/(lambda1 + x[i][1])
    ax1, ax2 = x[i][0]/(lambda2 + x[i][0]), x[i][1]/(lambda2 + x[i][1])
    gx2 = x[i][1]/(lambda3 + x[i][1])
    f1 = x[i][0]*(bx1 - ax1*x[i][0] - fx1*x[i][2])
    f2 = x[i][1]*(bx2 - ax2*x[i][1] - fx2*x[i][2])
    f3 = -x[i][2]*(mu + gx2*x[i][0] - fx1*x[i][1])

    x[i + 1][0] = x[i][0] + h*f1
    x[i + 1][1] = x[i][1] + h*f2
    x[i + 1][2] = x[i][2] + h*f3

x1, x2, x3 = [i[0] for i in x], [i[1] for i in x], [i[2] for i in x]

xc = [[x1c] * len(M), [x2c] * len(M), [x3c] * len(M)]
plt.figure(figsize=(10, 7))
plt.plot(M, x1, 'g', linewidth=3, label=r'$x_{1}$')
plt.plot(M, x2, 'b', linewidth=3, label=r'$x_{2}$')
plt.plot(M, x3, 'r', linewidth=3, label=r'$x_{3}$')
plt.plot(M, xc[0], 'k--', linewidth=3, label=r'$x_{1}^*$')
plt.plot(M, xc[1], 'k--', linewidth=3, label=r'$x_{2}^*$')
plt.plot(M, xc[2], 'k--', linewidth=3, label=r'$x_{3}^*$')
plt.style.use('seaborn-whitegrid')
plt.xlim(0, N)
plt.ylim(bottom=0)
plt.legend(loc='upper right')
plt.xlabel('Время, дни')
plt.ylabel('Популяция, ед/л')
plt.annotate(
    text,
    xy=(0, 1),
    textcoords='axes fraction',
    xytext=(0.01, 0.99),
    horizontalalignment='left',
    verticalalignment='top'
)

plt.show()