import numpy as np
import matplotlib.pyplot as plt
import os

def TimoshenkoBeamElement2D(Le, L_beam, E, A, I, phi):
    phi = phi*(L_beam/Le)**2 # change denominator to Le**2

    # Stiffness matrix
    K_beam = E*I/(Le**3*(1+phi)) * np.array([[0, 0, 0, 0, 0, 0], 
                                             [0, 12, 6*Le, 0, -12, 6*Le], 
                                             [0, 6*Le, (4+phi)*Le**2, 0, -6*Le, (2-phi)*Le**2], 
                                             [0, 0, 0, 0, 0, 0], 
                                             [0, -12, -6*Le, 0, 12, -6*Le], 
                                             [0, 6*Le, (2-phi)*Le**2, 0, -6*Le, (4+phi)*Le**2]])
    
    K_bar = E*A/Le * np.array([[1, 0, 0, -1, 0, 0], 
                               [0, 0, 0, 0, 0, 0], 
                               [0, 0, 0, 0, 0, 0], 
                               [-1, 0, 0, 1, 0, 0], 
                               [0, 0, 0, 0, 0, 0], 
                               [0, 0, 0, 0, 0, 0]])
    
    return K_beam + K_bar

def eliminate_dof(K, F, dof):
    K_reduced = np.delete(K, dof, axis=0)
    K_reduced = np.delete(K_reduced, dof, axis=1)
    F_reduced = np.delete(F, dof)

    return K_reduced, F_reduced

def w_analytical(x, F, L, E, I, phi):
    # If F > 0, the beam bends downwards. This is different convention from the implemented FEM
    # When phi=0, we recover the Euler-Bernoulli solution. 

    # Heaviside function: 1 for x >= L, else 0
    H = np.where(x >= L, 1.0, 0.0)

    ramp = (x - L) * H          # (x-L) * H(x-L)
    ramp3 = (x - L)**3 * H      # (x-L)^3 * H(x-L)

    bend = 1/6 * ramp3 + x**3/6 - L*x**2/2
    shear = -phi*L**2/12 * (ramp + x)

    w = F/(E*I) * (bend + shear)

    return w

if __name__ == "__main__":
    E = 70e9 # Pa
    L = 1.0 # m
    A = 0.01 # m^2
    I = 1e-6 # m^4

    N = 3 # number of elements
    F_mag = 1000 # tip load magnitude in N (downwards)

    plt.figure(figsize=(10, 6))
    line_color = ['#1f77b4', 'orange', 'red']
    k = 0
    for phi in [1.0, 3.0, 5.0]: # shear correction factor
        dx = L/N
        x = np.linspace(0, L, N+1) # node locations
        numDOF = 3*(N+1) # Each node has 3 DOF: x-disp, y-disp, rotation

        # Assemble global stiffness matrix
        K = np.zeros((numDOF, numDOF))
        for i in range(N):
            Ke_i = TimoshenkoBeamElement2D(dx, L, E, A, I, phi)
            n = Ke_i.shape[0] # Number of DOFs per element (6 for 2-node beam element)
            K[3*i:3*i+n, 3*i:3*i+n] += Ke_i

        # Load vector (point-load at tip)
        F = np.zeros(numDOF) #
        F[-2] = -F_mag # N
        
        # Apply boundary conditions
        dof_fixed = [0, 1, 2] # Fix the left end
        K_reduced, F_reduced = eliminate_dof(K, F, dof_fixed)

        # Solve for displacements
        u_reduced = np.linalg.solve(K_reduced, F_reduced)
        u = np.zeros(numDOF)
        u[3:] = u_reduced

        # Plot results
        ## Analytical solution for comparison
        x_ana = np.linspace(0, L, 200)
        w_ana = w_analytical(x_ana, F=F_mag, L=L, E=E, I=I, phi=phi)
        plt.plot(x_ana/L, w_ana/L, 'k-', linewidth=2, label='Analytical' if k == 0 else '_nolegend_', zorder=0)
        # FEM solution (y-displacement at nodes)
        w_nodes = u[1::3]
        plt.plot(x/L, w_nodes/L, '--o', color=line_color[k], zorder=1, label=r'FEM $\phi={}$'.format(phi))
        k += 1

    ## Euler Bernoulli for comparison
    x_EB = np.linspace(0, L, 200)
    w_EB = w_analytical(x_EB, F=F_mag, L=L, E=E, I=I, phi=0)
    plt.plot(x_EB/L, w_EB/L, 'g-', linewidth=2, label='Euler-Bernoulli', zorder=0)
    plt.xlabel(r'$x/L$', fontsize=14)
    plt.ylabel(r'$\frac{w(x)}{L}$', rotation=0, labelpad=10, fontsize=14)
    plt.suptitle('Cantilever Beam with Tip Load', fontsize=18)
    plt.title('Timoshenko vs. Euler-Bernoulli Beam (FEM N={})'.format(N), fontsize=16)
    plt.legend(fontsize=13)
    plt.grid()
    # plt.savefig('FEM_Solution.pdf')
    plt.savefig('FEM_Solution.png', dpi=400, bbox_inches='tight', pad_inches=0)

    # os.system('pdfcrop FEM_Solution.pdf')
    # os.system('mv FEM_Solution-crop.pdf FEM_Solution.pdf')