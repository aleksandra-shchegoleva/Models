import numpy as np
from typing import List
from math import exp, sin, pi, fabs
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.integrate import odeint
import seaborn as sns

mpl.rcParams.update({'font.size': 18})
np.random.seed(5)
u = [0]

class TPP_model_NAD():
    """TPP модель с непрерывным NAD-управлением"""
    __r, __K = 2, 108
    __a, __beta, __mu, __theta = 0.7, .6, .3654, 0.2
    __gamma = 5
    __eta, __k = 1, .1
    _xc, __T1 = (__gamma * __mu) / (__beta - __mu - __theta), 2
    __x0, N, __h = np.array([[50, 80]]), 100, 0.1
    __xi, __z0 = 0.2, 0
    __A, __w, __phi = 0, 0.1, 0
    __m_xi, __sko = 0, 1
    TEXT = ""

    def __init__(self, x0=[50, 80], z0=0, n=100, h=0.1):
        """Инициализация параметров модели"""
        self.__x0, self.N, self.__h = np.array([x0]), n, h
        self._m = np.arange(0, self.N, self.__h)
        self.__z0 = z0
        self.set_text()

    def set_text(self):
        self.TEXT = f"$x_{1}(0)$ = {self.__x0[0, 0]}\n" \
                    f"$x_{2}(0)$ = {self.__x0[0, 1]}\n" \
                    f"$z(0)$ = {self.__z0}\n" \
                    f"$r$ = {self.__r}\n" \
                    f"$K$ = {self.__K}\n" \
                    f"$\\alpha$ = {self.__a}\n" \
                    f"$\\beta$ = {self.__beta}\n" \
                    f"$\mu$ = {self.__mu}\n" \
                    f"$\\theta$ = {self.__theta}\n" \
                    f"$\gamma$ = {self.__gamma}\n" \
                    f"$\eta$ = {self.__eta}\n" \
                    f"$k$ = {self.__k}\n" \
                    f"$x_1^*$ = {self._xc}\n" \
                    f"$T_1$ = {self.__T1}\n" \
                    f"$M[\\xi]$ = {self.__m_xi}\n" \
                    f"$\\sigma$ = {self.__sko}"
                    # "Параметры шума\n" \
                    # f"$A$ = {self.__A}\n" \
                    # f"$w$ = {self.__w}\n" \
                    # f"$\\phi$ = {self.__phi}"

    def set_xi(self, xi):
        """Присвоение значения xi"""
        self.__xi = xi
        self.set_text()

    def set_A(self, A):
        self.__A = A
        self.set_text()

    def set_w(self, w):
        self.__w = w
        self.set_text()

    def set_phi(self, phi):
        self.__phi = phi
        self.set_text()

    def calc_xi(self, t):
        # return self.__A * sin(self.__w * t + self.__phi)
        # return self.__xi
        return self.__sko * np.random.randn()

    def set_sko(self, sko):
        self.__sko = sko
        self.set_text()

    def plot_solution(self, x, z, u):
        xc = np.array(self._xc * np.ones(len(self._m)))

        plt.figure(figsize=(10, 7))
        plt.plot(self._m, x[0], 'g', label=r'$x_{1}$')
        plt.plot(self._m, x[1], 'b', label=r'$x_{2}$')
        plt.plot(self._m, xc, 'k--', label=r'$x_{1}^*$')
        sns.set_style('whitegrid')
        plt.xlim(0, self.N)
        plt.ylim(bottom=0)
        plt.legend(loc="upper right")
        plt.xlabel('Время, дни')
        plt.ylabel('Популяция, ед/л')
        plt.annotate(
            self.TEXT,
            xy=(0, 1),
            textcoords='axes fraction',
            xytext=(0.99, 0.01),
            horizontalalignment='right',
            verticalalignment='bottom',
            fontsize=10
        )
        plt.figure(figsize=(10, 7))
        plt.plot(self._m, u, 'g', label=r'$u(t)$')
        sns.set_style('whitegrid')
        plt.xlim(0, self.N)
        plt.legend(loc="upper right")
        plt.xlabel('Время, дни')
        plt.ylabel('Управление')

        plt.figure(figsize=(10, 7))
        plt.plot(self._m, z, 'g', label=r'$z(t)$')
        sns.set_style('whitegrid')
        plt.xlim(0, self.N)
        plt.legend(loc="upper right")
        plt.xlabel('Время, дни')
        plt.ylabel('Возмущение')

        plt.show()

    def calc_euler_add(self, plot) -> List:
        x1, x2, z = [self.__x0[0, 0]].copy(), [self.__x0[0, 1]].copy(), [self.__z0].copy()
        for i in range(len(self._m) - 1):
            f = g = x1[i] / (self.__gamma + x1[i])

            psi = x1[i] - self._xc
            xi = self.calc_xi(self._m[i])
            z_dt = self.__eta * psi
            psi_end = psi + self.__k * z[i]
            u.append(- psi_end / self.__T1 - (self.__r * x1[i] * (1 - x1[i] / self.__K) - self.__a * f * x2[i]) - z[i] - self.__k * self.__eta * psi)
            f1 = self.__r * x1[i] * (1 - x1[i] / self.__K) - self.__a * f * x2[i] + u[-2]
            f2 = self.__beta * f * x2[i] - self.__mu * x2[i] - self.__theta * g * x2[i]

            x1.append(x1[-1] + self.__h * f1 + xi)
            x2.append(x2[-1] + self.__h * f2)
            z.append(z[-1] + self.__h * z_dt)

        if plot:
            self.plot_solution([x1, x2], z, u)

        return [self._m, [x1, x2], z]

    def get_eta(self):
        return self.__eta

    def get_k(self):
        return self.__k


def plot(sko):
    model.set_sko(sko)
    (M, x, z) = model.calc_euler_add(plot=True)
    eta, k = model.get_eta(), model.get_k()


model = TPP_model_NAD(n=100, h=0.1, z0=1, x0=[10, 80])
plot(sko=9.3)