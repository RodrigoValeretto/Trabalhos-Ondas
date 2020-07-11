import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

# Declaração de ctes
c = sp.speed_of_light
sigma = 0
mi = sp.mu_0  # 4*np.pi*(10**(-7))
epsilon = sp.epsilon_0  # 1/(mi*(c**2))
N = 150
I = 100
J = 100
d = 1
dt = 1.0005*d/(c*np.sqrt(2))
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


def unit_step(n):
    if n <= 20:
        return 1
    else:
        return 0


def YEE_2D(Ez, Hx, Hy):
    for n in range(1, N):
        for i in range(1, I):
            for j in range(1, J):
                Hx[n, i, j] = Da*Hx[n-1, i, j] + Db * \
                    (Ez[n-1, i, j] - Ez[n-1, i, j+1] - Mx*d)
                
                Hy[n, i, j] = Da*Hy[n-1, i, j] + Db * \
                    (Ez[n-1, i+1, j] - Ez[n-1, i, j] - My*d)

                Ez[n, i, j] = Ca*Ez[n-1, i, j] + Cb * \
                    (Hy[n, i, j] - Hy[n, i-1, j] +
                     Hx[n, i, j-1] - Hx[n, i, j] - Jz*d)

        Ez[n, int(I/2), int(J/2)] = unit_step(n)

    return Ez, Hx, Hy


Ez = np.zeros((N+1, I+1, J+1))      # Saltos n = 1/2, i = 1/2, j = 1/2
Hx = np.zeros((N+1, I+1, J+1))      # Saltos n = 1, i = 1/2, j = 1
Hy = np.zeros((N+1, I+1, J+1))      # Saltos n = 1, i = 1, j = 1/2

Ez, Hx, Hy = YEE_2D(Ez, Hx, Hy)

# Plot
fps = 30  # frame per sec
frn = N  # frame number of the animation
x = np.linspace(0, d*I, I)
y = np.linspace(0, d*J, J)
x, y = np.meshgrid(x, y)
zarray = np.zeros((N, I, J))

for n in range(frn):
    for i in range(I):
        for j in range(J):
            # Mude Ez para exibir o vetor desejado
            zarray[n, i, j] = Ez[n, i, j]


def update_plot(frame_number, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(
        x, y, zarray[frame_number, :, :], rstride=3, cstride=3, cmap="magma")


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plot = [ax.plot_surface(x, y, zarray[0, :, :], rstride=3, cstride=3)]
ax.set_zlim(-1, 1)

ani = animation.FuncAnimation(
    fig, update_plot, frn, fargs=(zarray, plot), interval=500/fps)

plt.show()
