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