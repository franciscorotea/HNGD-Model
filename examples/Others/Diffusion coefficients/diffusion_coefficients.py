import matplotlib.pyplot as plt
import numpy as np

'''A script to plot the different diffusion coefficients of Hydrogen in Zr
found in the literature.

References:
    
    [1] E.A. Gulbransen, K.F. Andrew, J. Electrochem. Soc. 101 (1954) 560–566.
    [2] M. Someno, Nippon Kinzoku Gakkaishi 24 (1960).
    [3] A. Sawatzky, J. Nucl. Mater. 2 (1960) 321–328.
    [5] J.J. Kearns, J. Nucl. Mater. 43 (1972) 330–338.
    [6] B. Kammenzind, D.G. Franklin, H.R. Peters, and W.J. Duffin, ASTM 
    Special Technical Publication 1295 (1996) 338–369.
    [7] D. Khatamian, J. Alloys Compd. 253–254 (1997) 471–474.
    [8] Y. Zhang, C. Jiang, X. Bai, Nature Sci. Rep. 7 (2017) 41033.
    
    '''

database = {'Gulbransen et al. (1954)': {'Reference': '[1]',
                                         'Material': 'Zirconium',
                                         'Isotope': 'Hydrogen',
                                         'D0': 1.09e-7,
                                         'Ed': 0.49},
            'Someno (1960)': {'Reference': '[2]',
                              'Material': 'Zirconium',
                              'Isotope': 'Hydrogen',
                              'D0': 4.15e-7,
                              'Ed': 0.41},
            'Sawatzky (1960)': {'Reference': '[3]',
                                'Material': 'Zircaloy-2',
                                'Isotope': 'Hydrogen',
                                'D0': 2.17e-7,
                                'Ed': 0.36},
            'Kearns (1972)': {'Reference': '[4]',
                              'Material': 'Zirconium/Zircaloy-2/Zircaloy-4',
                              'Isotope': 'Hydrogen',
                              'D0': 7.73e-7,
                              'Ed': 0.47},
            'Kammenzind et al. (1996)': {'Reference': '[5]',
                                         'Material': 'Zircaloy-4',
                                         'Isotope': 'Hydrogen',
                                         'D0': 0.8e-7,
                                         'Ed': 0.34},
            'Khatamian (1997)': {'Reference': '[6]',
                                 'Material': 'Zr-2.5%Nb',
                                 'Isotope': 'Deuterium',
                                 'D0': 2.61e-7,
                                 'Ed': 0.40},
            'Zhang et al. (2017)': {'Reference': '[7]',
                                    'Material': 'Zirconium',
                                    'Isotope': 'Hydrogen',
                                    'D0': 1.08e-6,
                                    'Ed': 0.46}
        }

k_B = 8.617333262e-5 # Boltzmann constant [eV K-1]

temperature = np.array([300, 1000]) + 273 # Temperature [K]

fig, ax = plt.subplots(figsize=(7,4))

fig.canvas.manager.set_window_title('Diffusion coefficient according to different authors')

for author, data in database.items():
    D = data['D0'] * np.exp(-data['Ed']/(k_B*temperature))
    ax.semilogy(10**4/temperature, D, label=author)
    
ax.legend(ncols=2, prop={'size': 8})

ax.set_xlabel('10$^4$/Temperature [K]')
ax.set_ylabel('Diffusion coefficient [m$^2$/s]')

sec_ax = ax.secondary_xaxis('top', functions=(lambda x: 1/(10**4/x), lambda x: 1/(10**4/x)))
sec_ax.set_xlabel('Temperature [°C]')
sec_ax.set_xticklabels(np.rint(10**4/ax.get_xticks() - 273).astype(int))

ax.set_ylim([1e-12, 1e-7])
ax.grid(which='both')

fig.tight_layout()

plt.savefig('diffusion_coefficients.png', bbox_inches='tight')

plt.show()
