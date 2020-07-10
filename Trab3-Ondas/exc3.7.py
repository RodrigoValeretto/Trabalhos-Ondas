import numpy as np
import scipy.constants as sp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.animation import FuncAnimation

# Declaração de ctes
c = sp.speed_of_light
sigma = 0
mi = 4*np.pi*(10**(-7))
epsilon = 1/(mi*(c**2))
N = 50
I = 50
J = 50
d = 0.05
dt = d/(c*np.sqrt(2))
X = d*I
Y = d*J
T = dt*N

# Cálculo das ctes
Jz = 0
Mx = 0
My = 0
Ca = (1 - sigma*dt/(2*epsilon))/(1 + sigma*dt/(2*epsilon))
Cb = (dt/(epsilon*d))/(1 + sigma*dt/(2*epsilon))
Da = (1 - sigma*dt/(2*mi))/(1 + sigma*dt/(2*mi))
Db = (dt/(mi*d))/(1 + sigma*dt/(2*mi))


def animate(n):
    line.set_data(Ez[n, x, y])
    return line,


def unit_step(t):
    if t >= 0:
        return 1
    else:
        return 0


def YEE_2D(Ez, Hx, Hy):
    for n in range(0, N):
        Ez[n, int(I/2), int(J/2)] = unit_step(n)

    for n in range(0, N):
        for i in range(0, I):
            for j in range(0, J-1):
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

x = np.linspace(0, d*I, I)
y = np.linspace(0, d*J, J)

fig = plt.figure()
ax = plt.axes(projection='3d')
line = ax.plot3D(x, y, Ez[0, np.arange(0, I), np.arange(0, J)])

# anim = FuncAnimation(fig, func=animate, frames=np.arange(0, N), interval=100, blit=True)

plt.show()
