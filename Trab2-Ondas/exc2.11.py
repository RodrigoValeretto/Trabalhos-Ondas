# IMPORTS
import numpy as np
import matplotlib.pyplot as plt

# Definição das constantes do exc
S1 = 1                     # Constante de courant para o primeiro plot
S2 = 1.075
K = 220                    # Número de divisões do espaço livre
N = K                      # Número de divisões do tempo
Nd1 = 190                  # N desejado 1
Nd2 = 200                  # N desejado 2


def func(n, S):
    if n > 74:
        return 0
    else:
        return np.exp(-((((n-60)/(5))**2)))


def calcWavef(U, S1, S2, N):
    # Definição dos valores iniciais
    for n in range(0, N):
        if n == 90:
            U[n, 0] = func(n, S2)
        else:
            U[n, 0] = func(n, S1)

    for i in range(0, K):
        U[0, i] = 0

    for n in range(1, N):
        for i in range(1, K):
            if i == 90:
                U[n+1, i] = (S2**2)*(U[n, i+1] - 2*U[n, i] +
                                     U[n, i-1]) + 2*U[n, i] - U[n-1, i]
            else:
                U[n+1, i] = (S1**2)*(U[n, i+1] - 2*U[n, i] +
                                     U[n, i-1]) + 2*U[n, i] - U[n-1, i]
    return U


# SCRIPT PRINCIPAL
U = np.zeros((N+1, K+1))        # Função de onda U1 no tempo n, no espaço i
gI = np.arange(0, 160, 1)       # Definição do grid i
gI2 = np.arange(70, 110, 1)     # Definição do grid com apenas 20 celulas
fig, ax = plt.subplots(2, 2)
U = calcWavef(U, S1, S2, N)

ax[0, 0].plot(gI, U[Nd1, gI])
ax[0, 0].grid(True)
ax[0, 0].set_title("Grid entre 0 e 220 e n = " + str(Nd1))
ax[0, 0].set_ylabel("Wavefunction U(i)")
ax[0, 0].set_ylim(-1, 1)

ax[0, 1].plot(gI, U[Nd2, gI], color='orange')
ax[0, 1].grid(True)
ax[0, 1].set_title("Grid entre 0 e 220 e n = " + str(Nd2))
ax[0, 1].set_ylabel("Wavefunction U(i)")
ax[0, 1].set_ylim(-1, 1)

ax[1, 0].plot(gI2, U[Nd1, gI2])
ax[1, 0].grid(True)
ax[1, 0].set_title("Grid entre 70 e 110 e n = " + str(Nd1))
ax[1, 0].set_ylabel("Wavefunction U(i)")
ax[1, 0].set_xlabel("Coordenada i do grid")
ax[1, 0].set_ylim(-1, 1)

ax[1, 1].plot(gI2, U[Nd2, gI2], color='orange')
ax[1, 1].grid(True)
ax[1, 1].set_title("Grid entre 70 e 110 e n = " + str(Nd2))
ax[1, 1].set_ylabel("Wavefunction U(i)")
ax[1, 1].set_xlabel("Coordenada i do grid")
ax[1, 1].set_ylim(-1, 1)

plt.show()
