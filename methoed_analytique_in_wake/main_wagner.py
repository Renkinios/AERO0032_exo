# Wagner theory of undsteady lift is develope :
# 1. In the time domaine 
# 2. On the baisis of the unsteady potential flow therory and the confoirmal Joukowski transformation 
# 3. The idea is to reprode tge effect of the starting vortex 


import numpy as np
import scalar_implement as si
import implement_matrix as im
import plot_fct as pf

c                      = 0.25  # [m]
thinkness              = 0.02  # [m]
rho_air                = 1.225 # [kg/m^3] at sea lvl
rho_alum               = 2700  # [kg/m^3]
mass                   = c * thinkness * rho_alum # [kg/m]
xf                     = c * 0.4 # [m]
omega_h_0              = 2 * np.pi * 2 # [rad/s]
omega_alpha_0          = 2 * np.pi * 8 # [rad/s]
psi_one                = 0.165
psi_two                = 0.335
eps_1                  = 0.0455
eps_2                  = 0.3


b             = c / 2
K_h           = si.get_structural_stifness_h(mass, omega_h_0)
inertia_alpha = si.get_inertia_alpha(mass, c, xf)
K_alpha       = si.get_structural_stifness_alpha(inertia_alpha, omega_alpha_0)
S             = si.get_moment_S(mass, c, xf)
e             = si.get_exentricite(c, xf)


def phi(t, U, b) : 
    return 1 - 0.165 * np.exp(-0.0455 * t * U/b) - 0.335 * np.exp(-0.3 * t *U/b)

def phi_dot(t, U, b, dt=1e-5):
    return (phi(t + dt, U, b) - phi(t, U, b)) / dt


U = np.linspace(5, 40, 100)

matrix_w = []
matrix_damping = []
for v in U : 
    M = im.get_M_wagner(mass, rho_air, b, S, xf, inertia_alpha, c)
    C = im.get_C_wagner(rho_air, phi(0, v, b), c, e, xf, v)
    K = im.get_K_wagner(K_h, K_alpha, rho_air, v, phi_dot(0, v, b), phi(0, v, b), xf, e, c)
    W = im.get_W_wagner(psi_one, psi_two, eps_1, eps_2, e, c, b, v, rho_air)
    W0 = im.get_W0_wagner(eps_1, eps_2, v, b)
    Q = im.get_Q_wagner(M, C, K, W, W0)
    eigenvalues = np.linalg.eigvals(Q)
    # print(eigenvalues)
    row_w = []
    row_damping = []
    for i in eigenvalues:
        row_w.append(np.abs(i))
        row_damping.append(-np.real(i)/np.abs(i))
    matrix_w.append(row_w)
    # print(matrix_w)
    matrix_damping.append(row_damping)
matrix_w = np.array(matrix_w).T
# print(matrix_w.size)
matrix_damping = np.array(matrix_damping).T

# pf.plot_wagner_w(matrix_w, U)
pf.plot_wagner_damping(matrix_damping, U)

    





