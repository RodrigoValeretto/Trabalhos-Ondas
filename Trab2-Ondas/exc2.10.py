# IMPORTS
import numpy as np
import matplotlib.pyplot as plt

# Definição das constantes do exc
S = 1.0005                 # Constante de courant para o primeiro plot
K = 220                    # Número de divisões do espaço livre
N = K                     # Número de divisões do tempo
Nd1 = 200                  # N desejado 1
Nd2 = 210                  # N desejado 2
Nd3 = 220                  # N desejado 3


def func(n, S):
    if n > 220/S:
        return 0
    else:
        return np.exp(-((((n-50/S)/(20/S))**2)))


def calcWavef(U, S, N):
    # Definição dos valores iniciais
    for n in range(0, N):
        U[n, 0] = func(n, S)

    for i in range(0, K):
        U[0, i] = 0

    for n in range(1, N):
        for i in range(1, K):
            U[n+1, i] = (S**2)*(U[n, i+1] - 2*U[n, i] +
                                U[n, i-1]) + 2*U[n, i] - U[n-1, i]
    return U


# SCRIPT PRINCIPAL
U = np.zeros((N+1, K+1))      # Função de onda U1 no tempo n, no espaço i
gI = np.arange(0, K, 1)     # Definição do grid i
fig, ax = plt.subplots(1, 3, sharey=True)
U = calcWavef(U, S, N)

ax[0].plot(gI, U[Nd1, gI])
ax[0].grid(True)
ax[0].set_title("Wavefunction para S = " + str(S) + " e n = " + str(Nd1))
ax[0].set_ylabel("Wavefunction U(i)")
ax[0].set_xlabel("Coordenada i do grid")

ax[1].plot(gI, U[Nd2, gI], color='orange')
ax[1].grid(True)
ax[1].set_title("Wavefunction para S = " + str(S) + " e n = " + str(Nd2))
ax[1].set_ylabel("Wavefunction U(i)")
ax[1].set_xlabel("Coordenada i do grid")

ax[2].plot(gI, U[Nd3, gI], color='purple')
ax[2].grid(True)
ax[2].set_title("Wavefunction para S = " + str(S) + " e n = " + str(Nd3))
ax[2].set_ylabel("Wavefunction U(i)")
ax[2].set_xlabel("Coordenada i do grid")

plt.show()
