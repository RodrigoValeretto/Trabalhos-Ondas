#imports
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
#import matplotlib

'''DEFINIÇÃO DAS CONSTANTES DO CÓDIGO'''
c = 3*(10**8)           #Valor constante que representa a velocidade da luz no vácuo
uf = 0.9*c              #Velocidade do sinal de tensão ou corrente
lmax = 30000             #Valor constante que define o tamanho da malha
test = (10+0.5)*lmax/uf       #Valor constante que define o momento de regime estacionário
tmax = test             #Valor constante que define o tempo máximo na malha
t1 = lmax/uf            #Tempo de trânsito, em que as ondas alcançam a carga
G = 0                   #Valor da condutância para linha sem perdas
R = 0                   #Valor da resistencia para linha sem perdas
Zo = 50                 #Valor da impedância para z=-lmax
L = Zo/uf               #Valor da indutância para linha sem perdas
C = 1/(Zo*uf)           #Valor da capacitância para linha sem perdas
Rs = 75                 #Valor da resistência da fonte especificada
Rl1 = np.inf            #Valor da carga no caso 1
Rl2 = 0                 #Valor da carga no caso 2
Rl3 = 100               #Valor da carga no caso 3
k = 100

'''DEFINIÇÃO DE FUNÇÕES UTILIZADAS'''
def u(t):               #Definição da função degrau de t
    if t < 0:
        return 0
    else:
        return 1

def vs1(t):              #Definição da função vs1 de t
    return 2*u(t)

def vs2(t):              #Definição da função vs2 de t
    return u(t) - u(t-lmax/(10*uf))

def restr():
    dz = lmax/k
    dt = dz/uf
    dt = dt*0.5
    n = int(tmax/dt)
    return dz,dt,n

def animateV(i):
    line.set_ydata(V[i,:])
    return line,

def animateI(i):
    line2.set_ydata(I[i,:])
    return line2,

def main(vs, Rl):
    #Definição da corrente e tensão na malha
    I = np.zeros((n, k-1))
    V = np.zeros((n, k))
    Va = np.zeros((n, k))

    '''CALCULO DE CONSTANTES'''
    c1 = (2*dt)/(Rs*C*dz)
    if Rl == 0:
        c2 = np.inf
    elif Rl == np.inf:
        c2 = 0
    else:
        c2 = (2*dt)/(Rl*C*dz)
    c3 = (dt**2)/(L*C*(dz**2))

    '''CONDIÇÕES DE CONTORNO E INICIAIS'''
    V[0, 0] = (Zo/(Rs + Zo))*vs(0)
    I[0, 0] = vs(0)/(Zo + Rs)

    #Va[0, 0] = C*(dz/dt)*V[0, 0]
    Va[0, 0] = V[0, 0]

    '''INICIO DAS ITERAÇÕES'''
    for j in range(1, n):
        Va[j, 0] = (1 - c1)*Va[j-1, 0] - 2*I[j-1, 0] + (2/Rs)*vs((j-1)*dt)

        Va[j, k-1] = (1 - c2)*Va[j-1, k-1] + 2*I[j-1, k-2]

        for i in range(1, k-1):
            Va[j, i] = Va[j-1, i] - (I[j-1, i] - I[j-1, i-1])
        for i in range(0, k-1):
            I[j, i] = I[j-1, i] - c3*(Va[j, i+1] - Va[j, i])

    '''for j in range(0, n):
        for i in range(0, k):
            V[j, i] = (dt/(C*dz))*Va[j, i]'''

    print("***Print da matriz de corrente***")
    print(I)
    print("*********************************")
    print("***Print da matriz de tensão***")
    print(Va)
    print("*******************************")

    '''    
    fig = plt.figure()
    x = np.arange(0, (k)*dz, dz)
    camera = Camera(fig)
    for i in range(0, n):
        plt.plot(x, Va[i])
        camera.snap()
    animation = camera.animate()
    plt.show()
    '''

    return Va, I

'''DEFINIÇÃO E PRINT DOS SALTOS, DIVISÕES E CTEs'''
dz, dt, n = restr()
print("***Condições e limites do programa***")
print("L =", L)
print("C =", C)
print("k =", k)
print("n =", n)
print("dz =", dz)
print("dt =", dt)
print("*************************************")

'''SCRIPT DO CÓDIGO PRINCIPAL'''
V, I = main(vs1, Rl3)
'''PLOT'''
x = np.linspace(0, lmax, k)
y = np.linspace(0, lmax, k-1)

fig, ax = plt.subplots()
line, = ax.plot(x, V[0])
plt.grid(True)
plt.xlim(-5, lmax+5)
plt.ylim(-0.2,0.2)
plt.title("Tensão")
anim = FuncAnimation(fig, func=animateV, frames=np.arange(0, n), interval=100, blit=True)

fig2, ax2 = plt.subplots()
line2, = ax2.plot(y, I[0])
plt.grid(True)
plt.xlim(-5, lmax+5)
plt.ylim(-0.2,0.2)
plt.title("Corrente")
anim = FuncAnimation(fig2, func=animateI, frames=np.arange(0, n), interval=100, blit=True)

plt.show()
#main(vs1, Rl2)
#main(vs1, Rl1)

