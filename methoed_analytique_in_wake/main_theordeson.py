# Implementation of the Theordeson function to take into acount the wake for aerolastic analysis
# The assumption made for this method : 
# 1. The flow is always attached.
# 2. The airfoil is a flat plat.
# 3. The wake is flat.

# Need to made the assumption that the flow is a sinusoidal flow this is true like this method going 
# to calculate juste the flutter speed, and at flutter is a sunoidal.
import numpy as np
import plot_fct as pf
import scalar_implement as si


c                      = 0.25  # [m]
thinkness              = 0.02  # [m]
rho_air                = 1.225 # [kg/m^3] at sea lvl
rho_alum               = 2700  # [kg/m^3]
mass                   = c * thinkness * rho_alum # [kg/m]
xf                     = c * 0.4 # [m]
omega_h_0              = 2 * np.pi * 2 # [rad/s]
omega_alpha_0          = 2 * np.pi * 8 # [rad/s]


b       = c / 2
K_h     = si.get_structural_stifness_h(mass, omega_h_0)
inertia = si.get_inertia_alpha(mass, c, xf)
K_alpha = si.get_structural_stifness_alpha(inertia, omega_alpha_0)
S       = si.get_moment_S(mass, c, xf)
e       = si.get_exentricite(c, xf)

# w = np.linspace(0, 100, 1000)
# U = np.linspace(0, 100, 1000)
# k = w * b / U

k = np.linspace(0.0001,10,100) # for the plot
Ck = 1 - (0.165)/(1-(0.0455/k) * 1j) - 0.335/(1 - (0.3/k) * 1j)

pf.plt_C_k(Ck, k)










