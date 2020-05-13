import math

'''DEFINIÇÃO DAS CONSTANTES DO CÓDIGO'''
c = 3*(10**8)           #Valor constante que representa a velocidade da luz no vácuo
lmax = 10               #Valor constante que define o tamanho da malha
tmax = 10               #Valor constante que define o tempo máximo na malha
test = 10*lmax/(0.9*c)  #Valor constante que define o momento de regime estacionário
n = 10                  #Valor constante que define o número de divisões na malha temporal(t)
k = 10                  #Valor constante que define o número de divisões na malha espacial(z)
L = 10                  #Valor da indutância fornecida pelo problema
C = 10                  #Valor da capacitância fornecida pelo problema

'''DEFINIÇÃO DE FUNÇÕES UTILIZADAS'''
def u(t):               #Definição da função degrau de t
    if t < 0:
        return 0
    else:
        return 1

def vs(t):              #Definição da função vs de t
    return 2*u(t)

def v(t):               #Definição da função v de t
    return u(t) - u(t-lmax/(10*0.9*c))

def telegraf(v, R, L):
    return 0

def restr(dz, dt, L, C):
    vp = 1/math.sqrt(L*C)
    if dt <= dz/vp:
        return True
    else:
        return False

'''SCRIPT DO CÓDIGO PRINCIPAL'''
dz = lmax/k
dt = tmax/n     #Pode ser reajustado com tempo estacionário futuramente

if not restr(dz, dt, L, C):
    print("Valores inválidos estipulados, modelo diverge")