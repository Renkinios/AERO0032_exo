def get_inertia_alpha(masse, chord, center_of_torsion) : 
    return 1/3 * masse * ( chord**2 - 3 * center_of_torsion * chord + 3 * center_of_torsion**2)

def get_moment_torsion(chord, mass, center_of_torsion) :
    return mass * (chord/2 - center_of_torsion)

def get_exentricite(chord, center_of_torsion) :
    return center_of_torsion/chord - 1/4

def get_structural_stifness_h(mass, omega_h) : 
    return omega_h**2 * mass

def get_structural_stifness_alpha(inertia, omega_alpha) : 
    return omega_alpha**2 * inertia