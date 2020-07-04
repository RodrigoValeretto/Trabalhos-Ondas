# IMPORTS
import numpy as np
import scipy.constants as sp

# SCRIPT
c = sp.speed_of_light   # Constante de velocidade da Luz
N = 40                  # Número de divisões do tempo
K = 40                  # Número de divisões do espaço livre
S = 1
lmax = 10               # Tamanho da malha de espaço (m)
tmax = 10               # Tamanho da malha de tempo (m)
dx = lmax/K             # Cálculo do valor do salto em espaço
dt = dx/(S*c)           # Cálculo do valor do salto em tempo
U = np.zeros((N, K))    # Função de onda U no tempo n, no espaço i

def func(t)
    return 1;

def calcWavef(U, S):
    # Definição dos valores iniciais
    for n in range(0, K):
        U(n, 0) = func(n*dt);

    if S == 1:          # Caso do salto de tempo mágico
        for n in range(0, N):
            for i in range(1, K):
                U(n+1, i) = U(n, i+1) + U(n, i-1) + U(n-1, i)
    else:
        for n in range(0, N):
            for i in range(1, K):
                U(n+1, i) = ((c*dt)**2)*((U(n, i+1) - 2*U(n, i) + U(n, i+1))/(dx**2)) + 2*U(n, i) - U(n-1, i)
    return U
