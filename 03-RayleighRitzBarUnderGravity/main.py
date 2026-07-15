import numpy as np
import matplotlib.pyplot as plt

def K_ab(a,b):
    p = a+b-1
    Kab = a*b*E*A*(L**p/p)
    return Kab

def F_a(a):
    Fa = rho*g*A*L**(a+1)/(a+1)
    return Fa

def u_h(u, x):
    n = len(u)
    u_h = 0
    for i in range(n):
        u_h += u[i] * x**(i+1)
    return u_h

if __name__ == "__main__":
    rho = 2700 # kg/m^3
    g = 9.81 # m/s^2
    E = 70e9 # Pa
    L = 1.0 # m
    A = 0.01 # m^2

    plt.figure(figsize=(8,6))
    for n in [1, 2, 3]: # order of the polynomial basis
        # Placeholder for the stiffness matrix and force vector
        K = np.zeros((n,n))
        F = np.zeros(n)

        # Construct the stiffness matrix and force vector
        for a in range(n):
            for b in range(n):
                K[a,b] = K_ab(a+1,b+1) # need +1 because equations were derived for a,b = 1,2,...,n
            F[a] = F_a(a+1)

        # Solve for the coefficients
        u = np.linalg.solve(K, F)
        print('Coefficients for n={}: {}'.format(n, u))

        # Choose a set of x values to evaluate the solution
        x = np.linspace(0, L, 100)

        # Plot the solution for this n
        ls = ':' if n==3 else '--'
        plt.plot(x, u_h(u, x), label='n={}'.format(n), linewidth=2, linestyle=ls, zorder=1)
        plt.title('Rayleigh-Ritz Bar Under Gravity (n={})'.format(n))
        plt.xlabel('Position along the bar (m)')
        plt.ylabel('Deflection (m)')
        plt.grid()

    # Plot the exact solution for comparison
    u_exact = rho*g/(2*E) * (2*L-x)*x
    plt.plot(x, u_exact, 'k-', label='Exact Solution', linewidth=3, zorder=0)
    plt.legend()

    # Save the figure
    plt.savefig('deflection.pdf'.format(n))