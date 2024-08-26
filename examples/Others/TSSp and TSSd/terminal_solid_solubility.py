import matplotlib.pyplot as plt
import numpy as np
from table_legend import tablelegend

'''

A script to plot the different terminal solid solubilities for dissolution
and precipitation found in the literature.

References:
    
    [1] Singh, R.N., Mukherjee, S., Gupta, A., et al., Terminal solid solubility 
    of hydrogen in Zr-alloy pressure tube material. J. Alloys Comp. 389, 
    102–112 (2005).
    [2] Giroldi, J.P., Vizcaíno, P., Flores, A.V., et al., Hydrogen terminal 
    solid solubility determinations in Zr–2.5Nb pressure tube microstructure in 
    an extended concentration range. J. Alloys Comp. 474, 140–146 (2009).
    [3] Une, K., Ishimoto, S., Terminal solid solubility of hydrogen in unalloyed 
    zirconium by differential scanning calorimetry. J. Nucl. Sci. Technol. 41, 
    949–952 (2004).
    [4] Tang, R., Yang, X., Dissolution and precipitation behaviors of hydrides 
    in N18, Zry-4 and M5 alloys. Int. J. Hydrogen Energy 34, 7269–7274 (2009).
    [5] Zanellato, O. et al., Synchrotron diffraction study of dissolution and 
    precipitation kinetics of hydrides in Zircaloy-4, J. Nucl. Mater. 420 (1–3)
    537–547 (2012).
    [6] Colas, K. B., Fundamental Experiments On Hydride Reorientation in 
    Zircaloy, The Pennsylvania State University (2012).
    [7] Slattery, G. F., The terminal solubility of hydrogen in zirconium alloys 
    between 30 and 400°C, J. Inst. Met. 95 (1967) 43.
    [8] McMinn, A., Darby, E. C., Schofield, J. S., Proceedings of the 12th 
    International Symposium on Zirconium in the Nuclear Industry, ASTM STP 1354, 
    pp.173–195 (2000).
    [9] Pan, Z. L., Ritchie, I. G., Puls, M. P., The terminal solid solubility 
    of hydrogen and deuterium in Zr-2.5Nb alloys, J. Nucl. Mater. 228, 227-237 
    (1996).
    [10] Khatamian, D., Pan, Z. L., Puls, M. P., Cann, C. D., Hydrogen solubility 
    limits in Excel, an experimental zirconium–based alloy, J. Alloys Compd. 231, 
    488–493 (1995).
    
'''

database = {'Singh et al. (Zry-2)': {'Reference': '[1]',
                                     'Year': '2005',
                                     'Material': 'Zircaloy-2',
                                     'Note': 'Cold-worked pressure tube',
                                     'TSSd0': 3.818e4,
                                     'TSSp0': 2.854e4,
                                     'Q_TSSd': 30004,
                                     'Q_TSSp': 25930,
                                     'linestyle': (0, (1, 1))},
# =============================================================================
#             'Singh et al. (Zr-2.5Nb)': {'Reference': '[1]',
#                                         'Year': '2005',
#                                         'Material': 'Zr-2.5Nb',
#                                         'Note': 'Cold-worked pressure tube',
#                                         'TSSd0': 6.358e4,
#                                         'TSSp0': 3.235e3,
#                                         'Q_TSSd': 35440,
#                                         'Q_TSSp': 17210,
#                                         'linestyle': (5, (10, 3))},
# =============================================================================
            'Giroldi et al. (Zr-2.5Nb)': {'Reference': '[2]',
                                          'Year': '2009',
                                          'Material': 'Zr-2.5Nb',
                                          'Note': 'Cold-worked pressure tube',
                                          'TSSd0': 5.26e4,
                                          'TSSp0': 2.22e4,
                                          'Q_TSSd': 31819,
                                          'Q_TSSp': 24618,
                                          'linestyle': (0, (3, 2))},
            'Une et al. (Zr)': {'Reference': '[3]',
                                'Year': '2004',
                                'Material': 'Zr',
                                'Note': '-',
                                'TSSd0': 1.41e5,
                                'TSSp0': 3.39e4,
                                'Q_TSSd': 38104,
                                'Q_TSSp': 27291,
                                'linestyle': (0, (3, 2, 1, 2))},
            'Tang et al. (N18)': {'Reference': '[4]',
                                  'Year': '2009',
                                  'Material': 'N18',
                                  'Note': '-',
                                  'TSSd0': 5.364e4,
                                  'TSSp0': 2.973e4,
                                  'Q_TSSd': 31809,
                                  'Q_TSSp': 25690,
                                  'linestyle': (0, (3, 1, 1, 1, 1, 1))},
            'Tang et al. (Zry-4)': {'Reference': '[4]',
                                    'Year': '2009',
                                    'Material': 'Zircaloy-4',
                                    'Note': '-',
                                    'TSSd0': 5.258e4,
                                    'TSSp0': 4.014e4,
                                    'Q_TSSd': 33117,
                                    'Q_TSSp': 27336,
                                    'linestyle': (0, (2, 1))},
            'Tang et al. (M5)': {'Reference': '[4]',
                                 'Year': '2009',
                                 'Material': 'M5',
                                 'Note': '-',
                                 'TSSd0': 8.497e4,
                                 'TSSp0': 3.064e4,
                                 'Q_TSSd': 34187,
                                 'Q_TSSp': 26180,
                                 'linestyle': (0, (1, 1, 2, 1, 2, 1))},
            'Zanellato et al. (Zry-4)': {'Reference': '[5]',
                                         'Year': '2012',
                                         'Material': 'Zircaloy-4',
                                         'Note': '-',
                                         'TSSd0': 510800,
                                         'TSSp0': 66440,
                                         'Q_TSSd': 45610,
                                         'Q_TSSp': 29630,
                                         'linestyle': (0, (1, 5))},
            'Colas (Zry-2)': {'Reference': '[6]',
                              'Year': '2012',
                              'Material': 'Zircaloy-2',
                              'Note': '-',
                              'TSSd0': 108150,
                              'TSSp0': 13320,
                              'Q_TSSd': 36360,
                              'Q_TSSp': 21170,
                              'linestyle': (0, (1, 3))},
            'Slattery (Zr-2.5Nb)': {'Reference': '[7]',
                                    'Year': '1967',
                                    'Material': 'Zr-2.5Nb',
                                    'Note': '-',
                                    'TSSd0': 6.860e4,
                                    'TSSp0': 4.110e4,
                                    'Q_TSSd': 33570,
                                    'Q_TSSp': 28000,
                                    'linestyle': (0, (2, 1, 1, 1, 1, 1, 1, 1))},
            'McMinn et al. (Zry-2)': {'Reference': '[8]',
                                      'Year': '2000',
                                      'Material': 'Zircaloy-2',
                                      'Note': '-',
                                      'TSSd0': 10.64e4,
                                      'TSSp0': 13.87e4,
                                      'Q_TSSd': 34629,
                                      'Q_TSSp': 34467,
                                      'linestyle': (0, (2, 3))},
            'Pan et al. (Zr-2.5Nb)': {'Reference': '[9]',
                                      'Year': '1996',
                                      'Material': 'Zr-2.5Nb',
                                      'Note': '-',
                                      'TSSd0': 8.080e4,
                                      'TSSp0': 2.812e4, # avg between 2.473e4 and 3.150e4
                                      'Q_TSSd': 34520,
                                      'Q_TSSp': 26915, # avg between 25840 and 27990
                                      'linestyle': (0, (2, 4))}  
# =============================================================================
#             'Khatamian et al. (Excel)': {'Reference': '[10]',
#                                          'Year': '1995',
#                                          'Material': 'Excel',
#                                          'Note': '-',
#                                          'TSSd0': 1.09e5,
#                                          'TSSp0': 1.45e5,
#                                          'Q_TSSd': 32600,
#                                          'Q_TSSp': 29900,
#                                          'linestyle': (0, (3, 2))}       
# =============================================================================
        }

R = 8.31446261815324                            # Gas constant [J mol-1 K-1]
temperature = np.linspace(0, 500, 500) + 273    # Temperature [K]

fig, ax = plt.subplots(figsize=(7,4))

fig.canvas.manager.set_window_title('TSSp and TSSd according to different authors')

for author, data in database.items():
    TSSp = data['TSSp0'] * np.exp(-data['Q_TSSp']/(R*temperature))
    ax.plot(temperature-273, TSSp, linestyle=data['linestyle'], color='#1f77b4', label=f'{author}, TSSp')
    
for author, data in database.items():
    TSSd = data['TSSd0'] * np.exp(-data['Q_TSSd']/(R*temperature))
    ax.plot(temperature-273, TSSd, linestyle=data['linestyle'], color='#ff7f0e', label=f'{author}, TSSd')
    
ax.set_xlabel('Temperature [°C]')
ax.set_ylabel('Terminal solid solubility [wt.ppm]')

ax.set_xlim([0, 500])
ax.set_ylim([0, 600])

tablelegend(ax, ncol=2,
            row_labels=[author for author in database],
            col_labels=['TSSp', 'TSSd'], 
            title_label='Author',
            prop={'size': 9})

fig.tight_layout()

plt.savefig('terminal_solid_solubility.png', bbox_inches='tight')

plt.show()
