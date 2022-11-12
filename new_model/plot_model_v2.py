import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import *
from typing import List, Container
import scipy

mpl.rcParams.update({'font.size': 18})


class Model:
    """Класс Model для поиска численного решения модели"""
    __B1, __B2 = 0.02, 0.03
    __A1, __A2 = 0.6, 0.45
    __e1, __e2 = 0.7, 0.85
    __K, __r, __mu = 10, 0.65, .001
    __x, N, __h = np.array([[1, 1, 1]]), 100, 0.1
    _m = np.arange(0, N, __h)
    TEXT = ""

    def __init__(self, x=[1, 1, 1], n=100, h=0.1,
                 B1=__B1, B2=__B2, A1=__A1, A2=__A2,
                 e1=__e1, e2=__e2, K=__K, r=__r, mu=__mu):
        """Инициализация параметров модели"""
        __B1, __B2 = 0.02, 0.03
        __A1, __A2 = 0.6, 0.45
        __e1, __e2 = 0.7, 0.85
        __K, __r, __mu = 10, 0.65, .001
        __x, __n, __h = np.array([x]), n, h
        _m = np.arange(0, __n, __h)

        TEXT = f"$x_{1}(0)$ = {__x[0, 0]}\n" \
               f"$x_{2}(0)$ = {__x[0, 1]}\n" \
               f"$x_{3}(0)$ = {__x[0, 2]}\n" \
               f"$r$ = {__r}\n" \
               f"$K$ = {__K}\n" \
               f"$\mu$ = {__mu},\n" \
               f"$B_{1}$ = {__B1}\n" \
               f"$B_{2}$ = {__B2}\n" \
               f"$A_{1}$ = {__A1}\n" \
               f"$A_{2}$ = {__A2}\n" \
               f"$e_{1}$ = {__e1}\n" \
               f"$e_{2}$ = {__e2}\n" \
               f"$h$ = {__h}"

    def calc_euler(self) -> Container:
        """Численное решение модели с помощью метода Эйлера"""
        """TODO: не работает"""
        # трофические функции b и f
        b = lambda x: self.__r * (1 - x / self.__K)
        f = lambda x1, x2: 1 / (1 + self.__B1 * x1 + self.__B2 * x2)

        for i in range(len(self._m) - 1):
            self.__x = np.append(self.__x, [[0, 0, 0]], axis=0)

            f1 = self.__x[i, 0] * (b(self.__x[i, 0]) - self.__A1 * f(self.__x[i, 0], self.__x[i, 1]) * self.__x[i, 2])
            f2 = self.__x[i, 1] * (b(self.__x[i, 1]) - self.__A2 * f(self.__x[i, 0], self.__x[i, 1]) * self.__x[i, 2])
            f3 = self.__x[i, 2] * (self.__mu + self.__e1 * f(self.__x[i, 0], self.__x[i, 1]) * self.__x[i, 0] -
                                   self.__e2 * f(self.__x[i, 0], self.__x[i, 1]) * self.__x[i, 1])

            self.__x[i + 1, 0] = self.__x[i, 0] + self.__h * f1
            self.__x[i + 1, 1] = self.__x[i, 1] + self.__h * f2
            self.__x[i + 1, 2] = self.__x[i, 2] + self.__h * f3
            print(self.__h * f3, f3)

        return self._m, self.__x

    def calc_runge_kutta(self):
        """Численное решение с помощью метода Рунге-Кутта"""


# cond1 = mu*(A1*B2*K + A1 - A2*B2*K)/(K*(A1 - A2))
# cond2 = mu*(-A1*B1*K + A2*B1*K + A2)/(K*(A1 - A2))
# cond3 = e1*(A1*B2*K + A1 - A2*B2*K)/(-A1*B1*K + A2*B1*K + A2)
# cond4= B1*mu + B2*mu + e1 + mu/K
# print(f"e_{2} = {cond1}")
# print(f"e_{1} = {cond2}")
# print(f"e_{2} = {cond3}")
# print(f"e_{2} = {cond4}")
# print(f"A_{1} = {A2*K*(B2*mu - e2)/(B2*K*mu - K*e2 + mu)}")
# print(f"A_{1} = {A2*(B1*K*mu + K*e1 + mu)/(K*(B1*mu + e1))}")
# print(f"A_{1} = {A2*(B1*K*e2 + B2*K*e1 + e2)/(B1*K*e2 + B2*K*e1 + e1)}")
# print(f"A_{2} = {A1*(B2*K*mu - K*e2 + mu)/(K*(B2*mu - e2))}")
# print(f"A_{2} = {A1*K*(B1*mu + e1)/(B1*K*mu + K*e1 + mu)}")
# print(f"A_{2} = {A1*(B1*K*e2 + B2*K*e1 + e1)/(B1*K*e2 + B2*K*e1 + e2)}")
# if e2 > cond1 and e1 > cond2 and e2 > cond3 and e2 > cond4:
#     x1c = (-A1 * B2 * K * mu + A1 * K * e2 - A1 * mu + A2 * B2 * K * mu - A2 * K * e2) / (
#             A1 * B1 * mu + A1 * e1 + A2 * B2 * mu - A2 * e2)
#     x2c = (A1 * B1 * K * mu + A1 * K * e1 - A2 * B1 * K * mu - A2 * K * e1 - A2 * mu) / (
#                 A1 * B1 * mu + A1 * e1 + A2 * B2 * mu - A2 * e2)
#     x3c = r * (B1 * K * mu + B2 * K * mu + K * e1 - K * e2 + mu) * (
#                 A1 * B1 * K * e2 + A1 * B2 * K * e1 + A1 * e1 - A2 * B1 * K * e2 - A2 * B2 * K * e1 - A2 * e2) / (
#                       K * (A1 * B1 * mu + A1 * e1 + A2 * B2 * mu - A2 * e2) ** 2)
# else:
#     x1c = K
#     x2c = K
#     x3c = 0

model = Model()
(M, x) = model.calc_euler()
print(x)

x1, x2, x3= [i[0] for i in x], [i[1] for i in x], [i[2] for i in x]

# xc = [[x1c] * len(M), [x2c] * len(M), [x3c] * len(M)]
# print(f"x_{1} = {x1c}, x_{2} = {x2c},x_{3} = {x3c}")
plt.figure(figsize=(10, 7))
plt.plot(M, x1, 'g', label=r'$x_{1}$')
plt.plot(M, x2, 'b', label=r'$x_{2}$')
plt.plot(M, x3, 'r', label=r'$x_{3}$')
# plt.plot(M, xc[0], 'k--', label=r'$x_{1}^*$')
# plt.plot(M, xc[1], 'k--', label=r'$x_{2}^*$')
# plt.plot(M, xc[2], 'k--', label=r'$x_{3}^*$')
plt.style.use('seaborn-whitegrid')
plt.xlim(0, model.N)
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.xlabel('Время, дни')
plt.ylabel('Популяция, ед/л')
# plt.annotate(
#     text,
#     xy=(0, 1),
#     textcoords='axes fraction',
#     xytext=(0.99, 0.01),
#     horizontalalignment='right',
#     verticalalignment='bottom',
#     fontsize=10
# )
# print(f"Параметры: x1(0) = {x1[0]}, x2(0) = {x2[0]}, x3(0) = {x3[0]}, \n r = {r}, K = {K}, mu = {mu}, lambda1 = {lambda1}, lambda3 = {lambda3} \n h = {h}")
plt.show()