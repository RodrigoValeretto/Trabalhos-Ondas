c = 3*(10**8)           #Valor constante que representa a velocidade da luz no vácuo 

def u(t):               #Definição da função degrau de t
    if t < 0:
        return 0
    else:
        return 1

def vs(t):              #Definição da função vs de t
    return 2*u(t)

def v(t):               #Definição da função v de t
    return u(t) - u(t-l/(10*0.9*c))

def telegraf(v, R, L):
