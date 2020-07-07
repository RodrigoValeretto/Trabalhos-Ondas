# IMPORTS
import numpy as np
import matplotlib.pyplot as plt

# Definição das constantes do exc
S1 = 1                              # Constante de courant para o primeiro plot
S2 = 0.25                           # Constante de courant para o segundo plot
S3 = 0.5                            # Constante de courant para o terceiro plot
K = 200                             # Número de divisões do espaço livre
Nd = 190                            # N desejado
N = int(K/((140*S1 + 60*S2)/200))    # Número de divisões do tempo
Nd = int(Nd/((140*S1 + 60*S2)/200))  # N desejado para fazer o plot 1


def func(n, S):
    if n > 100/S:
        return 0
    else:
        return np.exp(-((((n-50/S)/(20/S))**2)))


def calcWavef(U, S1, S2, N):
    # Definição dos valores iniciais
    for n in range(0, N):
        U[n, 0] = func(n, S1)

    for i in range(0, K):
        U[0, i] = 0

    for n in range(1, N-1):
        for i in range(1, K):
            if i < 140:
                U[n+1, i] = (S1**2)*(U[n, i+1] - 2*U[n, i] +
                                     U[n, i-1]) + 2*U[n, i] - U[n-1, i]
            else:
                U[n+1, i] = (S2**2)*(U[n, i+1] - 2*U[n, i] +
                                     U[n, i-1]) + 2*U[n, i] - U[n-1, i]
    return U


# SCRIPT PRINCIPAL
U = np.zeros((N, K+1))     # Função de onda U1 no tempo n, no espaço i
gI = np.arange(0, K, 1)     # Definição do grid i

U = calcWavef(U, S1, S2, N)
plt.plot(gI, U[Nd, gI])
plt.grid(True)
plt.title("Wavefunction para S = " + str(S1) + " e " + str(S2))
plt.ylabel("Wavefunction U(i)")
plt.xlabel("Coordenada i do grid")

plt.show()
