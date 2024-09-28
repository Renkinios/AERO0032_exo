import numpy as np
import matrix_implement as mt


def ComputePulsation(eigenvalues):
    omega = np.abs(eigenvalues)
    return omega

def ComputeDampingRatio(eigenvalues):
    damping = []
    for i in range(len(eigenvalues)):
        damping.append(-eigenvalues[i].real/np.abs(eigenvalues[i]))
    return damping


def FindFlutterSpeed(A, B, C, D, E, F, start_speed, rho_air) :
    tol = 10**-4
    iter_max = 1000
    i = 0
    U = start_speed
    delta_U = 1
    while i < iter_max :
        Q = mt.get_matrix_Q(A, B, C, D, E, F, U, rho_air)
        eigenvalues = np.linalg.eigvals(Q)
        zeta = ComputeDampingRatio(eigenvalues)
        # print("zeta : \t", zeta)
        if any(abs(z) < tol for z in zeta):
            return U
        
        elif any(z < 0 for z in zeta) : 
            delta_U = -delta_U/5
            U += delta_U
        else:
            U += delta_U
        
        
        
        i += 1
        
    return U
