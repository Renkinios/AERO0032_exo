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
import implement_matrix as im
import sympy as sp
from scipy.optimize import fsolve



c                      = 0.25  # [m]
thinkness              = 0.02  # [m]
rho_air                = 1.225 # [kg/m^3] at sea lvl
rho_alum               = 2700  # [kg/m^3]
mass                   = c * thinkness * rho_alum # [kg/m]
xf                     = c * 0.4 # [m]
omega_h_0              = 2 * np.pi * 2 # [rad/s]
omega_alpha_0          = 2 * np.pi * 8 # [rad/s]


b             = c / 2
K_h           = si.get_structural_stifness_h(mass, omega_h_0)
inertia_alpha = si.get_inertia_alpha(mass, c, xf)
K_alpha       = si.get_structural_stifness_alpha(inertia_alpha, omega_alpha_0)
S             = si.get_moment_S(mass, c, xf)
e             = si.get_exentricite(c, xf)

k = np.linspace(0.0001,10,100) # for the plot
Ck = 1 - (0.165)/(1-(0.0455/k) * 1j) - 0.335/(1 - (0.3/k) * 1j)

pf.plt_C_k(Ck, k)

w, U = sp.symbols('w U', real=True)
k = w * b / U
Ck = 1 - (0.165)/(1-(0.0455/k) * 1j) - 0.335/(1 - (0.3/k) * 1j)

D = im.Flutter_determinant(K_h, w, mass, rho_air, U, c, Ck, b, S, xf, K_alpha, inertia_alpha, e)
D = sp.Matrix(D)
print(D)

det_D = D.det()

# Séparer la partie réelle et imaginaire
real_part = sp.re(det_D)
imag_part = sp.im(det_D)

# 3. Transformer ces expressions symboliques en fonctions utilisables par fsolve
real_func = sp.lambdify((w, U), real_part, 'numpy')
imag_func = sp.lambdify((w, U), imag_part, 'numpy')

# 4. Définir la fonction d'équations pour fsolve
def equations(vars):
    w_val, U_val = vars
    return [real_func(w_val, U_val), imag_func(w_val, U_val)]

# 5. Résoudre avec une estimation initiale
initial_guess = [5, 25]  # Tu peux ajuster l'estimation initiale selon tes besoins
solution = fsolve(equations, initial_guess)

# 6. Afficher la solution
print(f'Solution: w = {solution[0]}, U = {solution[1]}')












