import matplotlib.pyplot as plt
import numpy as np
import matrix_implement as mt

def plot_pulsationAndDamping2speed(pitch_start, plunge_start, A, B, C, D, E, F, rho_speed, speed_flutter) :
    speed = np.linspace(5, 40, 100)
    pulsation_alpha = []
    pulsation_alpha.append(pitch_start)
    pulsation_h     = []
    pulsation_h.append(plunge_start)
    damping_alpha   = [0]
    damping_h       = [0]
    
    for V in speed :
        Q           = mt.get_matrix_Q(A, B, C, D, E, F, V, rho_speed)
        eigenvalues = np.linalg.eigvals(Q)
        omega_alpha = np.abs(eigenvalues[2])
        omega_h     = np.abs(eigenvalues[0])
        damp_alpha  = -eigenvalues[0].real/np.abs(eigenvalues[2])
        damp_h      = -eigenvalues[3].real/np.abs(eigenvalues[0])
        pulsation_alpha.append(omega_alpha) 
        pulsation_h.append(omega_h)
        damping_alpha.append(damp_alpha)
        damping_h.append(damp_h)
    speed = np.insert(speed, 0, 0)
    plt.figure(figsize=(8, 6))
    plt.plot(speed, pulsation_alpha, label="Pitch", color='b', linewidth=2)
    plt.plot(speed, pulsation_h, label="Plunge", color='r', linestyle='--', linewidth=2)
    plt.xlabel("Speed [m/s]", fontsize=12)
    plt.ylabel("Pulsation [rad/s]", fontsize=12)
    plt.title("Pulsation vs Speed", fontsize=14)
    plt.legend(loc="best")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("figure/pulsation.pdf", format="pdf", dpi=300)
    plt.close()

    # Plot for Damping Ratio
    plt.figure(figsize=(8, 6))
    plt.scatter(speed_flutter, 0, color='g', label="Flutter Speed")
    plt.plot(speed, damping_alpha, label="Pitch", color='b', linewidth=2)
    plt.plot(speed, damping_h, label="Plunge", color='r', linestyle='--', linewidth=2)
    plt.xlabel("Speed [m/s]", fontsize=12)
    plt.ylabel("Damping Ratio", fontsize=12)
    plt.title("Damping Ratio vs Speed", fontsize=14)
    plt.legend(loc="best")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("figure/damping.pdf", format="pdf", dpi=300)
    plt.close()
