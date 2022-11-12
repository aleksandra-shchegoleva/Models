import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import *

mpl.rcParams.update({'font.size': 18})

N = 500
h = 0.1
M = np.arange(0, N, h)
x = [[1, 1, 1]]
h = 0.1

B1, B2, A1, A2, e1, e2, K, r, mu = 0.02, 0.03, 0.6, 0.45, 0.7, 0.85, 10, 0.65, 0.001
cond1 = mu*(A1*B2*K + A1 - A2*B2*K)/(K*(A1 - A2))
cond2 = mu*(-A1*B1*K + A2*B1*K + A2)/(K*(A1 - A2))
cond3 = e1*(A1*B2*K + A1 - A2*B2*K)/(-A1*B1*K + A2*B1*K + A2)
cond4= B1*mu + B2*mu + e1 + mu/K
print(f"e_{2} = {cond1}")
print(f"e_{1} = {cond2}")
print(f"e_{2} = {cond3}")
print(f"e_{2} = {cond4}")
print(f"A_{1} = {A2*K*(B2*mu - e2)/(B2*K*mu - K*e2 + mu)}")
print(f"A_{1} = {A2*(B1*K*mu + K*e1 + mu)/(K*(B1*mu + e1))}")
print(f"A_{1} = {A2*(B1*K*e2 + B2*K*e1 + e2)/(B1*K*e2 + B2*K*e1 + e1)}")
print(f"A_{2} = {A1*(B2*K*mu - K*e2 + mu)/(K*(B2*mu - e2))}")
print(f"A_{2} = {A1*K*(B1*mu + e1)/(B1*K*mu + K*e1 + mu)}")
print(f"A_{2} = {A1*(B1*K*e2 + B2*K*e1 + e1)/(B1*K*e2 + B2*K*e1 + e2)}")
if e2 > cond1 and e1 > cond2 and e2 > cond3 and e2 > cond4:
    x1c = (-A1 * B2 * K * mu + A1 * K * e2 - A1 * mu + A2 * B2 * K * mu - A2 * K * e2) / (
            A1 * B1 * mu + A1 * e1 + A2 * B2 * mu - A2 * e2)
    x2c = (A1 * B1 * K * mu + A1 * K * e1 - A2 * B1 * K * mu - A2 * K * e1 - A2 * mu) / (
                A1 * B1 * mu + A1 * e1 + A2 * B2 * mu - A2 * e2)
    x3c = r * (B1 * K * mu + B2 * K * mu + K * e1 - K * e2 + mu) * (
                A1 * B1 * K * e2 + A1 * B2 * K * e1 + A1 * e1 - A2 * B1 * K * e2 - A2 * B2 * K * e1 - A2 * e2) / (
                      K * (A1 * B1 * mu + A1 * e1 + A2 * B2 * mu - A2 * e2) ** 2)
else:
    x1c = K
    x2c = K
    x3c = 0

text = f"$x_{1}(0)$ = {x[0][0]}\n" \
       f"$x_{2}(0)$ = {x[0][1]}\n" \
       f"$x_{3}(0)$ = {x[0][2]}\n" \
       f"$r$ = {r}\n" \
       f"$K$ = {K}\n" \
       f"$\mu$ = {mu},\n" \
       f"$B_{1}$ = {B1}\n" \
       f"$B_{2}$ = {B2}\n" \
       f"$A_{1}$ = {A1}\n" \
       f"$A_{2}$ = {A2}\n" \
       f"$e_{1}$ = {e1}\n" \
       f"$e_{2}$ = {e2}\n" \
       f"$h$ = {h}"

for i in range(len(M) - 1):
    x.append([0, 0, 0])
    bx1, bx2 = r*(1 - x[i][0]/K), r*(1 - x[i][1]/K)
    f = 1 / (1 + B1*x[i][0] + B2*x[i][1])
    f1 = x[i][0]*(bx1 - A1*f*x[i][2])
    f2 = x[i][1]*(bx2 - A2*f*x[i][2])
    f3 = -x[i][2]*(mu + e1*f*x[i][0] - e2*f*x[i][1])

    x[i + 1][0] = x[i][0] + h*f1
    x[i + 1][1] = x[i][1] + h*f2
    x[i + 1][2] = x[i][2] + h*f3

x1, x2, x3= [i[0] for i in x], [i[1] for i in x], [i[2] for i in x]

xc = [[x1c] * len(M), [x2c] * len(M), [x3c] * len(M)]
print(f"x_{1} = {x1c}, x_{2} = {x2c},x_{3} = {x3c}")
plt.figure(figsize=(10, 7))
plt.plot(M, x1, 'g', label=r'$x_{1}$')
plt.plot(M, x2, 'b', label=r'$x_{2}$')
plt.plot(M, x3, 'r', label=r'$x_{3}$')
plt.plot(M, xc[0], 'k--', label=r'$x_{1}^*$')
plt.plot(M, xc[1], 'k--', label=r'$x_{2}^*$')
plt.plot(M, xc[2], 'k--', label=r'$x_{3}^*$')
plt.style.use('seaborn-whitegrid')
plt.xlim(0, N)
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.xlabel('Время, дни')
plt.ylabel('Популяция, ед/л')
plt.annotate(
    text,
    xy=(0, 1),
    textcoords='axes fraction',
    xytext=(0.99, 0.01),
    horizontalalignment='right',
    verticalalignment='bottom',
    fontsize=10
)
# print(f"Параметры: x1(0) = {x1[0]}, x2(0) = {x2[0]}, x3(0) = {x3[0]}, \n r = {r}, K = {K}, mu = {mu}, lambda1 = {lambda1}, lambda3 = {lambda3} \n h = {h}")
plt.show()