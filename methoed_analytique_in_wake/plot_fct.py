import matplotlib.pyplot as plt
import numpy as np


def plt_C_k(C,k) :
    # Plot magnitude and phase of C(k)
    fig, axs = plt.subplots(2, 1, figsize=(8, 8))

    # Magnitude plot
    axs[0].plot(k, np.abs(C), color='b')
    axs[0].set_title('Magnitude of C(k)')
    axs[0].set_xlabel('k')
    axs[0].set_ylabel('|C(k)|')
    axs[0].grid(True)

    # Phase plot
    axs[1].plot(k, np.angle(C), color='r')
    axs[1].set_title('Phase of C(k)')
    axs[1].set_xlabel('k')
    axs[1].set_ylabel('Phase(C(k)) [radians]')
    axs[1].grid(True)

    plt.tight_layout()
    plt.savefig("figures/function_C.pdf", format='pdf', dpi=300)

def plot_wagner_w(W, U):
    fig, ax = plt.subplots()
    
    # Assurer que chaque courbe est tracée avec un label correspondant à son eigenvalue
    for idx, w in enumerate(W):
        ax.plot(U, w, label=f'Eigenvalue {idx+1}')
    
    ax.set_xlabel('Vitesse (m/s)')
    ax.set_ylabel('W (rad/s)')
    ax.grid(True)
    ax.legend()
    plt.show()

def plot_wagner_damping(damping, U):
    fig, ax = plt.subplots()
    
    # Assurer que chaque courbe est tracée avec un label correspondant à son eigenvalue
    for idx, d in enumerate(damping):
        ax.plot(U, d, label=f'Eigenvalue {idx+1}')
    # ax.plot(U, damping[0], label='Eigenvalue 1')
    # ax.plot(U, damping[1], label='Eigenvalue 2')
    # ax.plot(U, damping[2], label='Eigenvalue 1')
    # ax.plot(U, damping[3], label='Eigenvalue 2')
    # ax.plot(U, damping[4], label='Eigenvalue 3')
    # ax.plot(U, damping[5], label='Eigenvalue 4')
    # ax.plot(U, damping[6], label='Eigenvalue 1')
    # ax.plot(U, damping[7], label='Eigenvalue 2')
    # ax.plot(U, damping[8], label='Eigenvalue 3')


    
    ax.set_xlabel('Vitesse (m/s)')
    ax.set_ylabel('Damping ratio')
    ax.grid(True)
    ax.legend()
    plt.show()