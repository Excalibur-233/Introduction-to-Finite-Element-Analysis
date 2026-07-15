import numpy as np
import scipy.sparse as sp
from matplotlib import pyplot as plt

dx = 0.05
c = -3.8e-7 * dx**2
L = 1.0
N = int(L/dx) + 1
x = np.linspace(0, L, N)

k = [np.ones(N-1), -2*np.ones(N), np.ones(N-1)]
K = sp.diags(k, [-1, 0, 1], shape=(N, N))

K = K.toarray()

# enforces u(x=0) = 0
K[0,:] = np.zeros(N)
K[0,0] = 1

# enforces du/dx(x=L) = 0
K[-1,:] = np.zeros(N)
K[-1,-1] = 1
K[-1, -2] = -1

F = np.ones(N) * c
F[0] = 0 # Dirichlet BC at x=0
F[-1] = 0 # Neumann BC at x=L

u = np.linalg.solve(K, F)
u_analytical = c/(2*dx**2) * x * (x-2*L)

plt.figure()
plt.plot(x, u, 'r--', linewidth=2, label='Finite Difference')
plt.plot(x, u_analytical, 'k-', linewidth=2, label='Analytical')
plt.title("Displacement of bar under gravity")
plt.xlabel(r"$x$ [m]")
plt.ylabel("Displacement (m)")
plt.grid()
plt.legend()
plt.savefig("bar_under_gravity.pdf")