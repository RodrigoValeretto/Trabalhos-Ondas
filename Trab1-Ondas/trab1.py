#imports
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import random
#import matplotlib

'''DEFINIÇÃO DAS CONSTANTES DO CÓDIGO'''
c = 3*(10**8)           #Valor constante que representa a velocidade da luz no vácuo
uf = 0.9*c              #Velocidade do sinal de tensão ou corrente
lmax = 1000             #Valor constante que define o tamanho da malha
test = 10*lmax/uf       #Valor constante que define o momento de regime estacionário
tmax = test             #Valor constante que define o tempo máximo na malha
G = 0                   #Valor da condutância para linha sem perdas
R = 0                   #Valor da resistencia para linha sem perdas
Zo = 50                 #Valor da impedância para z=-lmax
L = Zo/uf               #Valor da indutância para linha sem perdas
C = 1/(Zo*uf)           #Valor da capacitância para linha sem perdas
Rs = 75                 #Valor da resistência da fonte especificada
Rl1 = math.inf          #Valor da carga no caso 1
Rl2 = 0                 #Valor da carga no caso 2
Rl3 = 100               #Valor da carga no caso 3

'''DEFINIÇÃO DE FUNÇÕES UTILIZADAS'''
def u(t):               #Definição da função degrau de t
    if t < 0:
        return 0
    else:
        return 1

def vs1(t):              #Definição da função vs1 de t
    return 2*u(t)

def vs2(t):              #Definição da função vs2 de t
    return u(t) - u(t-lmax/(10*0.9*c))

def restr(lmax, tmax, L, C):
    vp = 1/math.sqrt(L*C)
    k = random.randint(1, 200)
    n = random.randint(1, 200)
    dz = lmax/k
    dt = tmax/n
    while dt <= dz/vp:
        k = random.randint(1, 200)
        n = random.randint(1, 200)
        dz = lmax/k
        dt = tmax/n
    return dz,dt,k,n

'''SCRIPT DO CÓDIGO PRINCIPAL'''
dz, dt, k, n = restr(lmax, tmax, L, C)
print("***Condições e limites do programa***")
print("k =", k)
print("n =", n)
print("dz =", dz)
print("dt =", dt)
print("*************************************")

#Definição da corrente e tensão na malha
I = np.zeros((2*k, 2*n))
V = np.zeros((2*k, 2*n))

V[0][0:2*n] = Rl3/(Rl3 + Rs)*vs1(0)
I[0][0:2*n] = vs1(0)/(Rl3 + Rs)

#Definição dos c's usados nas fórmulas
c1 = -(2*dt)/(dt*dz*R + 2*dz*L)
c2 = (2*L - dt*R)/(2*L + dt*R)
c3 = -(2*dt)/(dt*dz*G + 2*dz*C)
c4 = (2*C - dt*G)/(2*C + dt*G)

#Definição das fórmulas pelo método FDTD
for j in range(0, 2*n-2):
    for i in range(0, 2*k-2):
        I[i][j+1] = c1*(V[i+1][j] - V[i-1][j]) + c2*I[i][j+1]
    for i in range(0, 2*k-2):
        V[i+1][j+2] = c3*(I[i+2][j+1] - I[i][j+1]) + c4*V[i+1][j]

print("***Print da matriz de corrente***")
print(I)
print("*********************************")
print("***Print da matriz de tensão***")
print(V)
print("*******************************")

'''TENTATIVA FALHA DE PLOT
x = np.arange(0, 2*k-2, 1)
y = np.arange(0, 2*n-2, 1)
z = I

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot_wireframe(x, y, I, color="b")
plt.show()
'''