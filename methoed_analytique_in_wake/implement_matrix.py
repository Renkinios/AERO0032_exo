import numpy as np

def Flutter_determinant(K_h, omega, mass, rho, U, c, Ck, b, S, xf, Kalpha, Ialpha, e) :
    M11 = (K_h - omega**2 * mass 
           +  np.pi * rho * U * c * Ck * omega *1j)
    
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