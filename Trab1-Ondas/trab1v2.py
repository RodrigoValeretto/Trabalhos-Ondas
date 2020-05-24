'''
Codigo referente ao trabalho 1 da disciplina de Ondas
Eletromagnéticas - SEL0612
Alunos responsáveis:
Rodrigo Augusto Valeretto       NUSP:10684792
Leonardo Cerce Guioto           NUSP:10716640
João Pedro Borges de Castro     NUSP:10276720
'''

#imports
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

'''DEFINIÇÃO DAS CONSTANTES DO CÓDIGO'''
c = 3*(10**8)           #Valor constante que representa a velocidade da luz no vácuo
uf = 0.9*c              #Velocidade do sinal de tensão ou corrente
lmax = 10000            #Valor constante que define o tamanho da malha
k = 100                 #Valor que representa o número de divisões da malha
test = (10)*lmax/uf     #Valor constante que define o momento de regime estacionário
tmax = test             #Valor constante que define o tempo máximo na malha
t1 = lmax/uf            #Tempo de trânsito, em que as ondas alcançam a carga
G = 0                   #Valor da condutância para linha sem perdas
R = 0                   #Valor da resistencia para linha sem perdas
Zo = 50                 #Valor da impedância para a linha
L = Zo/uf               #Valor da indutância para linha sem perdas
C = 1/(Zo*uf)           #Valor da capacitância para linha sem perdas
Rs = 75                 #Valor da resistência da fonte especificada
Rl1 = np.inf            #Valor da carga no caso 1
Rl2 = 0                 #Valor da carga no caso 2
Rl3 = 100               #Valor da carga no caso 3
gamac1 = (Rl1 - Zo)/(Rl1 + Zo)  #Coeficientes de reflexão das
gamac2 = (Rl2 - Zo)/(Rl2 + Zo)  #três diferentes cargas do
gamac3 = (Rl3 - Zo)/(Rl3 + Zo)  #problema, e também da fonte
gamag = (Rs - Zo)/(Rs + Zo) 


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

def restr():        #Função que calcula as restrições do problema
    dz = lmax/k     #A partir de lmax e k deduzimos todos os
    dt = dz/uf      #outros valores necessários
    dt = dt*0.5
    n = int(tmax/dt)
    return dz,dt,n

def animate(i): #Função responsável por gerar a animação do gráfico
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

def plot(i, j, ax, V, I, Vs, Rl):   #Função que faz o plot do gráfico
    linea, = ax[i,j].plot(x, V[0], color='b', label="Tensão")
    lineb, = ax[i,j].plot(y, I[0], color='y', label="Corrente")
    ax[i,j].set_title("Vs" + Vs + " e Rl=" + Rl)
    ax[i,j].set_xlim(-5, lmax+5)
    ax[i,j].set_ylim(fmin, fmax)
    ax[i,j].grid(True)
    ax[i,j].legend()
    return linea, lineb,


def main(vs, Rl):
    #Definição da corrente e tensão na malha
    I = np.zeros((n, k-1))
    V = np.zeros((n, k))
    Va = np.zeros((n, k))   #Tensão auxiliar

    '''CALCULO DE CONSTANTES'''
    c1 = (2*dt)/(Rs*C*dz)

    if Rl == 0:     #Condicional para as diversas resistencias
        c2 = np.inf
    elif Rl == np.inf:
        c2 = 0
    else:
        c2 = (2*dt)/(Rl*C*dz)

    c3 = (dt**2)/(L*C*(dz**2))

    '''CONDIÇÕES DE CONTORNO E INICIAIS'''
    V[0, 0] = (Zo/(Rs + Zo))*vs(0)
    I[0, 0] = vs(0)/(Zo + Rs)

    Va[0, 0] = C*(dz/dt)*V[0, 0]    #Multiplicação por constante para diminuir
                                    #os cálculos e poder aplicar o algoritmo

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
    '''FIM DAS ITERAÇÕES'''
    '''RECONVERSÃO DE VARIÁVEL AUXILIAR PARA TENSÃO ORIGINAL'''
    for j in range(0, n):
        for i in range(0, k):
            V[j, i] = (dt/(C*dz))*Va[j, i]

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
print("Gama da carga Rl1 =", gamac1)
print("Gama da carga Rl2 =", gamac2)
print("Gama da carga Rl3 =", gamac3)
print("Gama da Fonte =", gamag)
print("*************************************")

'''SCRIPT DO CÓDIGO PRINCIPAL'''
x = np.linspace(0, lmax, k)     #Essas variáveis são usadas em todos os plots
y = np.linspace(0, lmax, k-1)   #pois o dominio em z é igual para todos os casos
fmin = -3                       #Essas variáveis correspondem ao eixo vertical
fmax = 3                        #do grafico e podem ser alteradas para visualizar
                                #uma região especifica
fig, ax = plt.subplots(3,2, sharex='all')

'''CALCULANDO PARA VS1 E RL3'''
V1, I1 = main(vs1, Rl3)
'''PLOT VS1 E RL3'''
line, line2 = plot(0, 0, ax, V1, I1, "1", "100")

'''CALCULANDO PARA VS2 E RL3'''
V2, I2 = main(vs2, Rl3)
'''PLOT VS2 E RL3'''
line3, line4 = plot(0, 1, ax, V2, I2, "2", "100")

'''CALCULANDO PARA VS1 E RL2'''
V3, I3 = main(vs1, Rl2)
'''PLOT VS1 E RL2'''
line5, line6 = plot(1, 0, ax, V3, I3, "1", "0")

'''CALCULANDO PARA VS2 E RL2'''
V4, I4 = main(vs2, Rl2)
'''PLOT VS2 E RL2'''
line7, line8 = plot(1, 1, ax, V4, I4, "2", "0")

'''CALCULANDO PARA VS1 E RL1'''
V5, I5 = main(vs1, Rl1)
'''PLOT VS1 E RL1'''
line9, line10 = plot(2, 0, ax, V5, I5, "1", "inf")

'''CALCULANDO PARA VS2 E RL1'''
V6, I6 = main(vs2, Rl1)
'''PLOT VS2 E RL1'''
line11, line12 = plot(2, 1, ax, V6, I6, "2", "inf")

'''FUNÇÃO DE ANIMAÇÃO PARA TODOS OS PLOTS'''
anim = FuncAnimation(fig, func=animate, frames=np.arange(0, n), interval=100, blit=True)
#Caso queira AUMENTAR a velocidade de animação, basta DIMINUIR o valor de interval acima

'''FUNÇÃO QUE EXIBE O PLOT'''
plt.show()

'''FIM DO SCRIPT E DO CÓDIGO'''