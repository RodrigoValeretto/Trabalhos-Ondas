import scipy.constants as sp
import numpy as np
import matplotlib.pyplot as plt


#S = 1/np.sqrt(2)
S = 0.5
Nt = (2*np.pi*S)/(np.arccos(1-2*(S**2)))
c = sp.speed_of_light
lmax = 10
lmin = 1
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
        return 1
    else:
        return -zeta(N) - np.sqrt(zeta(N)**2 - 1)

'''Calculo dos valores da velocidade de fase e da cte de atenuação'''
N = 1
while N <= lmax:
    R1.append(calcVtil(N)/c)
    R2.append(calcAt(N))
    N += jump

'''PLOT DOS RESULTADOS'''
N = np.arange(lmin,lmax + jump, jump)

fig, ax = plt.subplots(1,2)
ax[0].plot(N, R1)
ax[1].plot(N, R2)
ax[0].grid(True)
ax[1].grid(True)
ax[0].set_title("Velocidade de fase normalizada")
ax[1].set_title("Cte de atenuação")
plt.show()

