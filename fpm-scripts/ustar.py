import numpy as np
import matplotlib.pyplot as plt


USTAR = np.array([4,5,6,7,8,9])
As_mennon = np.array([0.05, 0.58, 0.5, 0.39, 0.09, 0.05 ])
# As_fpm = np.array([0.0054, 0.545, 0.290, 0.140, 0.0642, 0.0631 ])
As_fpm = np.array([0.0054, 0.545, 0.410, 0.280, 0.0642, 0.0631 ])
St = 0.2 

fn = 0.2
fn = 0.75 * fn
# Ts = np.array([8.32, 6.78, 6.67, 6.23, 4.71, 4.85])
Ts = np.array([8.32, 6.78, 6.67, 6.0, 5.21, 4.55])
fs = 1 / Ts
fRatio = fs / fn

FIGSIZE_A = (5,3)
FIGSIZE_F = (5,2.2)
ANCHOR = (0,1.04)
PAD = 15
DPI = 500

SAVE_AMPLITUDE = "/home/simon/Desktop/VIV/Vstar-amplitude.png"
SAVE_FREQ = "/home/simon/Desktop/VIV/Vstar-freq.png"
plt.rcParams["mathtext.fontset"] = 'stix'

def main():

    plt.figure("Amplitude", figsize=FIGSIZE_A)

    plt.scatter(USTAR, As_mennon, marker='o', color = 'k', facecolors='none', label="Mennon & Mittal")
    plt.scatter(USTAR, As_fpm, marker='+', color='r', label="Present work")
    plt.plot([4], [0] , 'k--', label=r"Strouhal frequency $f_{St}$")
    plt.xlabel(rf"$V^*$")
    plt.ylabel(rf"$A_y^*$", rotation=0, labelpad=PAD)
    plt.legend(loc="lower left", frameon=False, bbox_to_anchor=ANCHOR)
    plt.tight_layout()
    plt.savefig(SAVE_AMPLITUDE, dpi=DPI)


    plt.figure("Frequency", figsize=FIGSIZE_F)
    plt.scatter(USTAR, fRatio, marker='+', color='r', label="Present work")
    plt.plot(USTAR, St * USTAR, 'k--', label=r"Strouhal frequency $f_{St}$")
    plt.xlabel(r"$V^*$")
    plt.ylabel(r"$\dfrac{f}{f_n}$", rotation=0, labelpad=PAD)
    # plt.legend(loc="lower left", frameon=False, bbox_to_anchor=ANCHOR)
    plt.tight_layout()
    plt.savefig(SAVE_FREQ, dpi=DPI)
    return


if __name__ == "__main__":
    main()

    plt.show()
