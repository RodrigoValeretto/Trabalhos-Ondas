#imports
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

'''DEFINIÇÃO DAS CONSTANTES DO CÓDIGO'''
c = 3*(10**8)           #Valor constante que representa a velocidade da luz no vácuo
uf = 0.9*c              #Velocidade do sinal de tensão ou corrente
lmax = 30000             #Valor constante que define o tamanho da malha
test = (10)*lmax/uf       #Valor constante que define o momento de regime estacionário
tmax = test             #Valor constante que define o tempo máximo na malha
t1 = lmax/uf            #Tempo de trânsito, em que as ondas alcançam a carga
G = 0                   #Valor da condutância para linha sem perdas
R = 0                   #Valor da resistencia para linha sem perdas
Zo = 50                 #Valor da impedância para z=-lmax
L = Zo/uf               #Valor da indutância para linha sem perdas
C = 1/(Zo*uf)           #Valor da capacitância para linha sem perdas
Rs = 75                 #Valor da resistência da fonte especificada
Rl1 = np.inf               #Valor da carga no caso 1
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

def animate(i):
    line.set_ydata(V1[i,:])
    line2.set_ydata(I1[i,:])
    line3.set_ydata(V2[i,:])
    line4.set_ydata(I2[i,:])
    line5.set_ydata(V3[i,:])
    line6.set_ydata(I3[i,:])
    line7.set_ydata(V4[i,:])
    line8.set_ydata(I4[i,:])
    line9.set_ydata(V5[i,:])
    line10.set_ydata(I5[i,:])
    line11.set_ydata(V6[i,:])
    line12.set_ydata(I6[i,:])
    return line, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12
'''
def animateI(i):
    line2.set_ydata(I[i,:])
    return line2,
'''
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

    Va[0, 0] = C*(dz/dt)*V[0, 0]
    #Va[0, 0] = V[0, 0]

    '''INICIO DAS ITERAÇÕES'''
    for j in range(1, n):
        Va[j, 0] = (1 - c1)*Va[j-1, 0] - 2*I[j-1, 0] + (2/Rs)*vs((j-1)*dt)

        if Rl == 0:         #Condição de curto-circuito
            Va[j, k-1] = 0
        else:               #Condição para outros valores
            Va[j, k-1] = (1 - c2)*Va[j-1, k-1] + 2*I[j-1, k-2]

        for i in range(1, k-1):
            Va[j, i] = Va[j-1, i] - (I[j-1, i] - I[j-1, i-1])
        for i in range(0, k-1):
            I[j, i] = I[j-1, i] - c3*(Va[j, i+1] - Va[j, i])

    for j in range(0, n):
        for i in range(0, k):
            V[j, i] = (dt/(C*dz))*Va[j, i]
    '''
    print("***Print da matriz de corrente***")
    print(I)
    print("*********************************")
    print("***Print da matriz de tensão***")
    print(Va)
    print("*******************************")
    
        
    fig = plt.figure()
    x = np.arange(0, (k)*dz, dz)
    camera = Camera(fig)
    for i in range(0, n):
        plt.plot(x, Va[i])
        camera.snap()
    animation = camera.animate()
    plt.show()
    '''

    return V, I

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
'''CALCULANDO PARA VS1 E RL3'''
V1, I1 = main(vs1, Rl3)

'''PLOT VS1 E RL3'''
x = np.linspace(0, lmax, k)     #Essas variáveis são usadas em todos os plots
y = np.linspace(0, lmax, k-1)   #pois o dominio em z é igual para todos os casos
fmin = -3
fmax = 3

fig, ax = plt.subplots(3,2)
line, = ax[0,0].plot(x, V1[0], color='b', label="Tensão")
line2, = ax[0,0].plot(y, I1[0], color='y', label="Corrente")
ax[0,0].set_title("Vs1 e Rl=100")
ax[0,0].set_xlim(-5, lmax+5)
ax[0,0].set_ylim(fmin, fmax)
ax[0,0].grid(True)
ax[0,0].legend()

'''CALCULANDO PARA VS2 E RL3'''
V2, I2 = main(vs2, Rl3)

'''PLOT VS2 E RL3'''
line3, = ax[0,1].plot(x, V2[0], color='b', label="Tensão")
line4, = ax[0,1].plot(y, I2[0], color='y', label="Corrente")
ax[0,1].set_title("Vs2 e Rl=100")
ax[0,1].set_xlim(-5, lmax+5)
ax[0,1].set_ylim(fmin, fmax)
ax[0,1].grid(True)
ax[0,1].legend()

'''CALCULANDO PARA VS1 E RL2'''
V3, I3 = main(vs1, Rl2)

'''PLOT VS1 E RL2'''
line5, = ax[1,0].plot(x, V3[0], color='b', label="Tensão")
line6, = ax[1,0].plot(y, I3[0], color='y', label="Corrente")
ax[1,0].set_title("Vs1 e Rl=0")
ax[1,0].set_xlim(-5, lmax+5)
ax[1,0].set_ylim(fmin, fmax)
ax[1,0].grid(True)
ax[1,0].legend()

'''CALCULANDO PARA VS2 E RL2'''
V4, I4 = main(vs2, Rl2)

'''PLOT VS2 E RL2'''
line7, = ax[1,1].plot(x, V4[0], color='b', label="Tensão")
line8, = ax[1,1].plot(y, I4[0], color='y', label="Corrente")
ax[1,1].set_title("Vs2 e Rl=0")
ax[1,1].set_xlim(-5, lmax+5)
ax[1,1].set_ylim(fmin, fmax)
ax[1,1].grid(True)
ax[1,1].legend()

'''CALCULANDO PARA VS1 E RL1'''
V5, I5 = main(vs1, Rl1)

'''PLOT VS1 E RL1'''
line9, = ax[2,0].plot(x, V5[0], color='b', label="Tensão")
line10, = ax[2,0].plot(y, I5[0], color='y', label="Corrente")
ax[2,0].set_title("Vs1 e Rl=inf")
ax[2,0].set_xlim(-5, lmax+5)
ax[2,0].set_ylim(fmin, fmax)
ax[2,0].grid(True)
ax[2,0].legend()

V6, I6 = main(vs2, Rl1)

'''PLOT VS1 E RL1'''
line11, = ax[2,1].plot(x, V6[0], color='b', label="Tensão")
line12, = ax[2,1].plot(y, I6[0], color='y', label="Corrente")
ax[2,1].set_title("Vs2 e Rl=inf")
ax[2,1].set_xlim(-5, lmax+5)
ax[2,1].set_ylim(fmin, fmax)
ax[2,1].grid(True)
ax[2,1].legend()

'''FUNÇÃO DE ANIMAÇÃO PARA TODOS OS PLOTS'''
anim = FuncAnimation(fig, func=animate, frames=np.arange(0, n), interval=100, blit=True)

plt.show()

