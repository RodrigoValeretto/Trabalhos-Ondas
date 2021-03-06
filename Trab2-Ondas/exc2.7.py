# IMPORTS
import numpy as np
import matplotlib.pyplot as plt

# Definição das constantes do exc
S1 = 1
S2 = 0.99
S3 = 0.5
K = 200                          # Número de divisões do espaço livre
N1 = int(K/S1)                 # Número de divisões do tempo
N2 = int(K/S2)                 # Número de divisões do tempo
N3 = int(K/S3)                 # Número de divisões do tempo
Nd1 = int(160/S1)                # N desejado para fazer o plot
Nd2 = int(160/S2)                # N desejado para fazer o plot
Nd3 = int(160/S3)                # N desejado para fazer o plot


def func(n, S):
    if n > 40/S:
        return 0
    else:
        return 1


def calcWavef(U, S, N):
    # Definição dos valores iniciais
    for n in range(0, N):
        U[n, 0] = func(n, S)

    for i in range(0, K):
        U[0, i] = 0

    for n in range(1, N-1):
        for i in range(1, K):
            U[n+1, i] = (S**2)*(U[n, i+1] - 2*U[n, i] +
                                U[n, i-1]) + 2*U[n, i] - U[n-1, i]
    return U


# SCRIPT PRINCIPAL
U1 = np.zeros((N1, K+1))    # Função de onda U1 no tempo n, no espaço i
U2 = np.zeros((N2, K+1))    # Função de onda U no tempo n, no espaço i
U3 = np.zeros((N3, K+1))    # Função de onda U no tempo n, no espaço i
gI = np.arange(0, K, 1)
fig, ax = plt.subplots(1, 3, sharey=True)

U1 = calcWavef(U1, S1, N1)
ax[0].plot(gI, U1[Nd1, gI])
ax[0].grid(True)
ax[0].set_title("Wavefunction para S = " + str(S1))
ax[0].set_ylabel("Wavefunction U(i)")
ax[0].set_xlabel("Coordenada i do grid")

U2 = calcWavef(U2, S2, N2)
ax[1].plot(gI, U2[Nd2, gI], color='orange')
ax[1].grid(True)
ax[1].set_title("Wavefunction para S = " + str(S2))
ax[1].set_ylabel("Wavefunction U(i)")
ax[1].set_xlabel("Coordenada i do grid")

U3 = calcWavef(U3, S3, N3)
ax[2].plot(gI, U3[Nd3, gI], color='purple')
ax[2].grid(True)
ax[2].set_title("Wavefunction para S = " + str(S3))
ax[2].set_ylabel("Wavefunction U(i)")
ax[2].set_xlabel("Coordenada i do grid")

plt.show()
