import numpy as np

def Flutter_determinant(K_h, omega, mass, rho, U, c, Ck, b, S, xf, Kalpha, Ialpha, e) :
    M11 = (K_h - omega**2 * mass 
           +  np.pi * rho * U * c * Ck * omega *1j
           - omega**2 * rho * np.pi * b**2)
    
    M12 = (-omega**2 * S 
       + rho * np.pi * b**2 * U * 1j * omega 
       + rho * np.pi * b**2 * (xf - c/2) * omega**2 
       + np.pi * rho * U * c * Ck * (U + (3/4 * c - xf) * omega * 1j))
    
    M21 = (-omega**2 * S 
           - np.pi * rho * U * e * c**2 * Ck * omega * 1j
           + (xf - c/2) * rho * np.pi * b**2 * omega**2)
    M22 = (Kalpha 
           - omega**2 * Ialpha 
           + (3/4 * c - xf) * rho * np.pi * b**2 * U * omega * 1j
           - np.pi * rho * U * e * c**2 * Ck * (U + (3/4 * c - xf) * omega * 1j)
           - (xf - c/2)**2 * rho * np.pi * b**2 * omega**2
           - (rho * np.pi * b**4)/8 * omega**2)
    return np.array([[M11, M12], [M21, M22]])

def get_M_wagner(m, rho, b, S, x_f, Ialpha, c):
    M11 = m + rho * np.pi * b**2
    M12 = S - rho * np.pi * b**2 * (x_f - c/2)
    M21 = S - rho * np.pi * b**2 * (x_f - c/2)
    M22 = Ialpha + rho * np.pi * b**2 *((x_f - c/2)**2 + b**2/8)
    return np.array([[M11, M12], [M21, M22]])

def get_C_wagner(rho, phi_0, c, e, x_f, U):
    C11 = phi_0
    C12 = c/4 + phi_0 * (3 * c/4 - x_f)
    C21 = - e * c * phi_0
    C22 = (3 * c/4 - x_f) * (c/4 - e * c * phi_0)
    return  np.pi * rho * U * c * np.array([[C11, C12], [C21, C22]])

def get_K_wagner(Kh, Kalpha ,rho, U, phi_dot_0, phi_0, x_f, e, c):
    K11 = Kh + np.pi * rho * U * c * phi_dot_0
    K12 = np.pi * rho * U * c * (U * phi_0 + (3 * c/4 - x_f) * phi_dot_0)
    K21 = - np.pi * rho * U * e * c**2 * phi_dot_0
    K22 = Kalpha - np.pi * rho * U * e * c**2 *(U * phi_0 + (3*c/4 - x_f) * phi_dot_0) 
    return np.array([[K11, K12], [K21, K22]])

def get_W_wagner(psi_one, psi_two, eps_1, eps_2, e, c, b, U, rho) :
    W11 = - psi_one * eps_1**2/b
    W12 = - psi_two * eps_2**2/b
    W13 = psi_one * eps_1 * (1 - eps_1 * (1 - 2 * e))
    W14 = psi_two * eps_2 * (1 - eps_2 * (1 - 2 * e))
    W21 = e * c * psi_one * eps_1**2/b
    W22 = e * c * psi_two * eps_2**2/b
    W23 = - e * c * psi_one * eps_1 * (1 - eps_1 * (1 - 2 * e))
    W24 = - e * c * psi_two * eps_2 * (1 - eps_2 * (1 - 2 * e))
    return  2 * np.pi * rho * U**3 * np.array([[W11, W12, W13, W14], [W21, W22, W23, W24]])


def get_W0_wagner(eps_1, eps_2, U, b) :
    W0_1 = [1, 0 , - eps_1 * U/b, 0, 0, 0]
    W0_2 = [1, 0 , 0, - eps_2 * U/b, 0, 0]
    W0_3 = [0, 1 , 0, 0, - eps_1 * U/b, 0]
    W0_4 = [0, 1 , 0, 0, 0, - eps_2 * U/b]
    return np.array([W0_1, W0_2, W0_3, W0_4])


def get_Q_wagner(M, C, K, W, W_0) :

    M_inv = np.linalg.inv(M)
    Q11   = - M_inv @ C
    Q12   = - M_inv @ K
    Q13   = - M_inv @ W    
    Q21   = np.eye(2)
    Q22   = np.zeros((2,2))
    Q23   = np.zeros((2,4))
    Q31   = np.zeros((4,2))
    Q32   = W_0
    return np.vstack([np.block([Q11, Q12, Q13]), np.block([Q21, Q22, Q23]), np.block([Q31, Q32])])