import scipy.constants as sp
import numpy as np
import matplotlib.pyplot as plt


S = 0.5
Nt = (2*np.pi*S)/(np.arccos(1-2*(S**2)))
c = sp.speed_of_light
lmax = 80
lmin = 3
jump = 0.1
R1 = []
R2 = []

'''FUNÇÕES UTILIZADAS'''


def calcVtil(N):
    if N >= Nt:
        return (2*np.pi*c)/(N*np.arccos(1 + 4*(np.cos(np.pi/N) - 1)))
    else:
        return (2*c)/N


def zeta(N):
    return 1 + ((1/S)**2)*(np.cos(2*np.pi*S/N) - 1)


def calcAt(N):
    if N >= Nt:
        return 0
    else:
        return -np.log(-zeta(N) - np.sqrt(zeta(N)**2 - 1))


'''Calculo dos valores da velocidade de fase e da cte de atenuação'''
N = lmin
while N <= lmax:
    R1.append((1-calcVtil(N)/c)*100)
    N += jump

'''PLOT DOS RESULTADOS'''
N = np.arange(lmin, lmax + jump, jump)

plt.plot(N, R1)
plt.grid(True)
plt.title("Porcentagem da Velocidade de fase relativa a velocidade da luz")
plt.show()
