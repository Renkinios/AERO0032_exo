
import matrix_implement as mt
import scalar_implement as st
import FindFlutterSpeed as ff
import numpy as np

c                      = 0.25  # [m]
thinkness              = 0.02  # [m]
rho_air                = 1.225 # [kg/m^3] at sea lvl
rho_alum               = 2700  # [kg/m^3]
mass                   = c * thinkness * rho_alum # [kg/m]
xf                     = c * 0.4 # [m]
omega_h_0              = 2 * np.pi * 2 # [rad/s]
omega_alpha_0          = 2 * np.pi * 8 # [rad/s]

# Scalar implementation
inertia_alpha             = st.get_inertia_alpha(mass, c, xf)
moment_torsion            = st.get_moment_torsion(c, mass, xf)
exentricity               = st.get_exentricite(c, xf)
structural_stifness_h     = st.get_structural_stifness_h(mass, omega_h_0)
structural_stifness_alpha = st.get_structural_stifness_alpha(inertia_alpha, omega_alpha_0)

matrix_A = mt.get_matrix_A(mass, c, inertia_alpha)
matrix_B = mt.get_matrix_B(c, xf)
matrix_D = mt.get_matrix_D(c, xf, exentricity)
matrix_E = mt.get_matrix_E(structural_stifness_h, structural_stifness_alpha)
matrix_F = mt.get_matrix_F(c, exentricity)
matrix_C = np.zeros((2, 2))


start_speed = 5 
speed_fluter = ff.FindFlutterSpeed(matrix_A, matrix_B, matrix_C, matrix_D, matrix_E, matrix_F, start_speed, rho_air)
print("speed fluter : \t", speed_fluter)


