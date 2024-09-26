import numpy as np

def get_matrix_A(m, S, inertia_alpha) : 
    return np.array([[m, S], [S, inertia_alpha]])

def get_matrix_B(chord, center_of_torsion) : 
    return np.pi * chord/2 * np.array([[1 , (chord/2 - center_of_torsion)], 
                                       [(chord/2 - center_of_torsion), 
                                        (chord/2 - center_of_torsion)**2 + (chord/2)**2 / 8]])

def get_matrix_D(chord, center_of_torsion, exentricity) :
    return chord * np.pi * np.array([[1, (3 * chord/4 - center_of_torsion) + chord/4], 
                                     [ - exentricity * chord, (chord/2 - center_of_torsion)**2 
                                      +  (3 * chord/4 - center_of_torsion) * chord/4]])

def get_matrix_E(structural_stifness_h, structural_stifness_alpha) : 
    return np.array([[structural_stifness_h, 0], [0, structural_stifness_alpha]])

def get_matrix_F(chord, exentricity) : 
    return chord * np.pi * np.array([[0, 1], [0, -exentricity * chord]])

def get_matrix_Q(A, B, C, D, E, F, speed, rho):
    M = A + rho * B
    inverse_M = np.linalg.inv(M)  # Inversion de la matrice A
    return np.block([[-inverse_M @ (C + rho * speed * D), -inverse_M @ (E + rho * speed**2 * F)], 
                     [np.eye(2), np.zeros((2, 2))]])

