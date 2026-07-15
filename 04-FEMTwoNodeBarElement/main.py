import numpy as np
import matplotlib.pyplot as plt
import os

def Ke(E,A,dx):
    return E*A/dx * np.array([[1, -1], [-1, 1]])

def eliminate_dof(K, F, dof):
    K_reduced = np.delete(K, dof, axis=0)
    K_reduced = np.delete(K_reduced, dof, axis=1)
    F_reduced = np.delete(F, dof)
    return K_reduced, F_reduced

def u_reconstruction(x, u, dx, numptsperelement=20):
    # Reconstruction of continuous u(x) for linear elements
    N = len(u) - 1 # Number of elements
    x_dense = []
    u_dense = []
    for i in range(N):
        x_i = x[i]
        x_ip1 = x[i + 1]
        u_i = u[i]
        u_ip1 = u[i + 1]

        # Linear interpolation
        x_local = np.linspace(0, dx, numptsperelement, endpoint=False)
        u_dense.append((1-x_local/dx)*u_i + (x_local/dx)*u_ip1)
        x_dense.append(x_i + x_local)
    
    # Add the last point
    x_dense.append([x[-1]])
    u_dense.append([u[-1]])

    return np.concatenate(u_dense), np.concatenate(x_dense)

if __name__ == "__main__":
    rho = 2700 # kg/m^3
    g = 9.81 # m/s^2
    E = 70e9 # Pa
    L = 1.0 # m
    A = 0.01 # m^2

    plt.figure(figsize=(10, 6))
    line_color = ['#1f77b4', 'orange', 'red']
    k = 0
    for N in [1, 2, 5]:
        dx = L/N
        x = np.linspace(0, L, N+1)

        # Assemble global stiffness matrix
        K = np.zeros((N+1, N+1))
        for i in range(N):
            Ke_i = Ke(E, A, dx)
            K[i:i+2, i:i+2] += Ke_i

        # Load vector (self-weight)
        w = rho * g * A # N/m
        F = np.zeros(N+1)
        for i in range(N):
            F[i] += w * dx / 2
            F[i+1] += w * dx / 2
        
        # Apply boundary conditions
        dof_fixed = [0] # Fix the left end
        K_reduced, F_reduced = eliminate_dof(K, F, dof_fixed)

        # Solve for displacements
        u_reduced = np.linalg.solve(K_reduced, F_reduced)
        u = np.zeros(N+1)
        u[1:] = u_reduced

        # Reconstruct continuous u(x) for plotting
        u_recon, x_recon = u_reconstruction(x, u, dx)

        # Plot results
        plt.plot(x_recon/L, u_recon/L, '-', linewidth=2, label='FEM N={}'.format(N), zorder=1, color=line_color[k])
        plt.scatter(x/L, u/L, marker='o', color=line_color[k], zorder=1)
        k += 1

    # Analytical solution for comparison
    x_analytical = np.linspace(0, L, 1000)
    u_analytical = -rho * g/(2*E) * x_analytical * (x_analytical - 2*L)
    plt.plot(x_analytical/L, u_analytical/L, 'k-', linewidth=3, label='Analytical Solution', zorder=0)
    plt.xlabel(r'$x/L$', fontsize=14)
    plt.ylabel(r'$\frac{u(x)}{L}$', rotation=0, labelpad=10, fontsize=14)
    plt.title('Bar Under Self-Weight', fontsize=18)
    plt.legend(fontsize=13)
    plt.grid()
    plt.savefig('FEM_Solution.pdf')

    os.system('pdfcrop FEM_Solution.pdf')
    os.system('mv FEM_Solution-crop.pdf FEM_Solution.pdf')