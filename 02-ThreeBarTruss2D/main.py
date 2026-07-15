import numpy as np
import matplotlib.pyplot as plt

def Rot(phi):
    R = np.array([[np.cos(phi), np.sin(phi), 0, 0],
                  [0, 0, np.cos(phi), np.sin(phi)]])
    return R

def tilde_K(E, A, L):
    # 1D stiffness matrix for a bar element
    K = E * A / L * np.array([[1, -1], [-1, 1]])
    return K

def K_e(R, tildeK):
    # General 2D stiffness matrix for a bar element
    return R.T @ tildeK @ R

def eliminate_DOFs(K,F,DOF_to_eliminate):
    # Eliminate specified DOFs from the stiffness matrix and load vector
    K_reduced = np.delete(K, DOF_to_eliminate, axis=0)
    K_reduced = np.delete(K_reduced, DOF_to_eliminate, axis=1)
    F_reduced = np.delete(F, DOF_to_eliminate)
    return K_reduced, F_reduced

if __name__ == "__main__":
    # Material and geometric properties
    E = 210e9  # Young's modulus in Pascals
    A = 0.01   # Cross-sectional area in square meters
    L = 1.0    # Length of each bar in meters

    # Angles of the bars with respect to the horizontal axis (in radians)
    phi1 = 0          # Bar 1 angle
    phi2 = 2*np.pi/3  # Bar 2 angle
    phi3 = np.pi/3    # Bar 3 angle

    # Rotation matrices for each bar element
    R1 = Rot(phi1)
    R2 = Rot(phi2)
    R3 = Rot(phi3)

    # Stiffness matrices for each bar element in axial coordinates
    tildeK1 = tilde_K(E, A, L)
    tildeK2 = tilde_K(E, A, L)
    tildeK3 = tilde_K(E, A, L)

    # Stiffness matrices for each bar element in 2D coordinates
    K1 = K_e(R1, tildeK1)
    K2 = K_e(R2, tildeK2)
    K3 = K_e(R3, tildeK3)

    # Stiffness matrices for each bar element in global coordinates
    K1_global = np.zeros((6, 6))
    K1_global[0:4, 0:4] = K1 # bar 1 connects node 1 (DOFs 0-1) and node 2 (DOFs 2-3)
    
    K2_global = np.zeros((6, 6))
    K2_global[2:6, 2:6] = K2 # bar 2 connects node 2 (DOFs 2-3) and node 3 (DOFs 4-5)
    
    K3_global = np.zeros((6, 6))
    K3_global[0:2, 0:2] = K3[0:2, 0:2] # bar 3 connects node 1 (DOFs 0-1) and node 3 (DOFs 4-5)
    K3_global[4:6, 4:6] = K3[2:4, 2:4]
    K3_global[0:2, 4:6] = K3[0:2, 2:4]
    K3_global[4:6, 0:2] = K3[2:4, 0:2]

    # Global stiffness matrix assembly
    K_global = K1_global + K2_global + K3_global

    # Define the load vector (assuming a load of -1000 N in the y-direction at node 3)
    F = np.zeros(6)
    F[5] = -1000  # Load at node 3 in the y-direction
    
    # Eliminate DOFs in global stiffness matrix
    DOF_to_eliminate = [1, 3, 4] # y-dir for node 1 and 2, x-dir for node 3
    K_reduced, F_reduced = eliminate_DOFs(K_global, F, DOF_to_eliminate=DOF_to_eliminate) # Eliminate DOFs for node 1 and node 2

    # Solve reduced system
    u_reduced = np.linalg.solve(K_reduced, F_reduced)

    U = np.zeros(6)
    U[DOF_to_eliminate] = 0
    U[[i for i in range(6) if i not in DOF_to_eliminate]] = u_reduced
    
    F_ext = K_global @ U
    print("F_ext = ", F_ext)

    ## Plot
    plt.figure(figsize=(8, 6))
    # Original positions of the nodes
    nodes = np.array([[0, 0], [L, 0], [L/2, L*np.sqrt(3)/2]])
    plt.plot(nodes[:, 0], nodes[:, 1], 'ko', label='Original Structure')
    plt.plot([nodes[0, 0], nodes[1, 0]], [nodes[0, 1], nodes[1, 1]], 'k-') # Bar 1
    plt.plot([nodes[1, 0], nodes[2, 0]], [nodes[1, 1], nodes[2, 1]], 'k-') # Bar 2
    plt.plot([nodes[0, 0], nodes[2, 0]], [nodes[0, 1], nodes[2, 1]], 'k-') # Bar 3
    # Deformed positions of the nodes
    scale = 5e5
    deformed_nodes = nodes + U.reshape(3, 2) * scale
    plt.plot(deformed_nodes[:, 0], deformed_nodes[:, 1], 'ro', label='Deformed Structure')
    plt.plot([deformed_nodes[0, 0], deformed_nodes[1, 0]], [deformed_nodes[0, 1], deformed_nodes[1, 1]], 'r-') # Bar 1
    plt.plot([deformed_nodes[1, 0], deformed_nodes[2, 0]], [deformed_nodes[1, 1], deformed_nodes[2, 1]], 'r-') # Bar 2
    plt.plot([deformed_nodes[0, 0], deformed_nodes[2, 0]], [deformed_nodes[0, 1], deformed_nodes[2, 1]], 'r-') # Bar 3
    plt.title('Three Bar Truss Deformation (Scale = {:e})'.format(scale))
    plt.xlabel('X-axis (m)')
    plt.ylabel('Y-axis (m)')
    plt.legend()
    plt.grid()
    plt.axis('equal')
    plt.savefig('symmetric_deformation.pdf')

    ## Alternatively, we can have a nonsymmetric deformation by not fixing u_1^3 to be 0
    # Let us set u_1^1 = 0 instead
    DOF_to_eliminate = [0, 1, 3] # x-dir for node 1, y-dir for node 1 and node 2
    K_reduced, F_reduced = eliminate_DOFs(K_global, F, DOF_to_eliminate=DOF_to_eliminate)
    u_reduced = np.linalg.solve(K_reduced, F_reduced)
    U = np.zeros(6)
    U[DOF_to_eliminate] = 0
    U[[i for i in range(6) if i not in DOF_to_eliminate]] = u_reduced
    
    F_ext = K_global @ U
    print("F_ext = ", F_ext)

    # Plot
    plt.figure(figsize=(8, 6))
    # Original positions of the nodes
    nodes = np.array([[0, 0], [L, 0], [L/2, L*np.sqrt(3)/2]])
    plt.plot(nodes[:, 0], nodes[:, 1], 'ko', label='Original Structure')
    plt.plot([nodes[0, 0], nodes[1, 0]], [nodes[0, 1], nodes[1, 1]], 'k-') # Bar 1
    plt.plot([nodes[1, 0], nodes[2, 0]], [nodes[1, 1], nodes[2, 1]], 'k-') # Bar 2
    plt.plot([nodes[0, 0], nodes[2, 0]], [nodes[0, 1], nodes[2, 1]], 'k-') # Bar 3
    # Deformed positions of the nodes
    scale = 5e5
    deformed_nodes = nodes + U.reshape(3, 2) * scale
    plt.plot(deformed_nodes[:, 0], deformed_nodes[:, 1], 'ro', label='Deformed Structure')
    plt.plot([deformed_nodes[0, 0], deformed_nodes[1, 0]], [deformed_nodes[0, 1], deformed_nodes[1, 1]], 'r-') # Bar 1
    plt.plot([deformed_nodes[1, 0], deformed_nodes[2, 0]], [deformed_nodes[1, 1], deformed_nodes[2, 1]], 'r-') # Bar 2
    plt.plot([deformed_nodes[0, 0], deformed_nodes[2, 0]], [deformed_nodes[0, 1], deformed_nodes[2, 1]], 'r-') # Bar 3
    plt.title('Three Bar Truss Deformation (Scale = {:e})'.format(scale))
    plt.xlabel('X-axis (m)')
    plt.ylabel('Y-axis (m)')
    plt.legend()
    plt.grid()
    plt.axis('equal')
    plt.savefig('nonsymmetric_deformation.pdf')