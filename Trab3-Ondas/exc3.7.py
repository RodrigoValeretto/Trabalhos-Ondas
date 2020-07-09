import numpy as np
import scipy.constants as sp

# Declaração de ctes
c = sp.speed_of_light
sigma = 1
mi = 4*np.pi*(10**(-7))
epsilon = 1/(mi*(c**2))
N = 100
I = 100
J = 100
d = 10
dt = d/(c*np.sqrt(2))

# Cálculo das ctes
Jz = 0
Mx = 0
My = 0
Ca = (1 - sigma*dt/(2*epsilon))/(1 + sigma*dt/(2*epsilon))
Cb = (dt/(epsilon*d))/(1 + sigma*dt/(2*epsilon))
Da = (1 - sigma*dt/(2*mi))/(1 + sigma*dt/(2*mi))
Db = (dt/(mi*d))/(1 + sigma*dt/(2*mi))


def YEE_2D(Ez, Hx, Hy):
    for n in range(0, N):
        for i in range(0, I):
            for j in range(0, J):
                Ez[n+1, i-1, j+1] = Ca*Ez[n-1, i-1, j+1] + Cb * \
                    (Hy[n, i, j+1] - Hy[n, i-1, j+1] +
                     Hx[n, i-1, j] - Hx[n, i-1, j+1] - Jz*d)

                Hx[n+1, i-1, j+1] = Da*Hx[n, i-1, j+1] + Db * \
                    (Ez[n+1, i-1, j+1] - Ez[n+1, i-1, j+3] - Mx*d)

                Hy[n+1, i, j+1] = Da*Hy[n, i, j+1] + Db * \
                    (Ez[n+1, i+1, j+1] - Ez[n+1, i-1, j+1] - My*d)
    return Ez, Hx, Hy


# Script Principal
Ez = np.zeros((N+1, I+1, J+3))      # Saltos n = 1/2, i = 1/2, j = 1/2
Hx = np.zeros((N+1, I+1, J+1))      # Saltos n = 1, i = 1/2, j = 1
Hy = np.zeros((N+1, I+1, J+1))      # Saltos n = 1, i = 1, j = 1/2

Ez, Hx, Hy = YEE_2D(Ez, Hx, Hy)