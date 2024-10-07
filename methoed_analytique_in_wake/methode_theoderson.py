import sympy as sp
import numpy as np
from scipy.optimize import fsolve

# 1. Définir les variables symboliques
x, y = sp.symbols('x y')

# 2. Définir une fonction complexe (par exemple det_D)
det_D = (x + sp.I*y)**2 + 1  # Juste un exemple de fonction complexe

# 3. Séparer la partie réelle et imaginaire
real_part = sp.re(det_D)
imag_part = sp.im(det_D)

# 4. Transformer les parties réelle et imaginaire en fonctions numériques
real_func = sp.lambdify((x, y), real_part, 'numpy')
imag_func = sp.lambdify((x, y), imag_part, 'numpy')

# 5. Définir une fonction qui combine les deux parties pour fsolve
def equations(vars):
    x_val, y_val = vars
    # Nous voulons que les deux parties soient nulles
    return [real_func(x_val, y_val), imag_func(x_val, y_val)]

# 6. Utiliser fsolve avec une approximation initiale
initial_guess = [10, 25]
solution = fsolve(equations, initial_guess)

print(f'Solution: x = {solution[0]}, y = {solution[1]}')
