import numpy as np
import matplotlib.pyplot as plt
import hngd.hngd
import hngd.default_model_parameters

"""

Reproduce Fig 4 from Passelaigue F., Simon, P. A. and Motta, A. T., 'Predicting 
the hydride rim by improving the solubility limits in the Hydride 
Nucleation-Growth-Dissolution (HNGD) model', Journal of Nuclear Materials 558, 
pp. 153363 (2022).
        
"""

model_parameters = hngd.default_model_parameters.model_parameters

TSSd0 = model_parameters['TSSd0']
Q_TSSd = model_parameters['Q_TSSd']
TSSp0 = model_parameters['TSSp0']
Q_TSSp = model_parameters['Q_TSSp']

temperature = 600
Cp = np.linspace(0, 16000, 1000)

TSSp, TSSd = hngd.hngd.calc_terminal_solid_solubility(TSSp0, Q_TSSp, TSSd0, Q_TSSd, temperature)
Cp_at = hngd.hngd.ppm_to_atomic_fraction(Cp)
C_delta = -9.93e-11*temperature**3 + 8.48e-8*temperature**2 - 5.73e-5*temperature + 0.623
TSSd_at = hngd.hngd.ppm_to_atomic_fraction(TSSd)
lvl_rule_Cp = Cp_at / (C_delta - TSSd_at)

def modifying_polinomial(g, delta):
    return g*lvl_rule_Cp - ((1-delta)*TSSd + g)*lvl_rule_Cp**2

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

for g in [100,150,200,250]:
    TSSd_eff = TSSd + modifying_polinomial(g=g, delta=1.05)
    ax1.plot(Cp, TSSd_eff, label=f'g = {g}')

ax1.set_xlabel('Hydride content Cp [wt.ppm]')
ax1.set_ylabel('Hydrogen effective solubility TSSd$_{eff}$ [wt.ppm]')
ax1.legend()

for delta in [0.85,0.95,1.05,1.15]:
    TSSd_eff = TSSd + modifying_polinomial(g=200, delta=delta)
    ax2.plot(Cp, TSSd_eff, label=f'delta = {delta}')

ax2.set_xlabel('Hydride content Cp [wt.ppm]')
ax2.set_ylabel('Hydrogen effective solubility TSSd$_{eff}$ [wt.ppm]')
ax2.legend()

fig.tight_layout()

plt.savefig('test_TSSd_eff_mHNGD.png', bbox_inches='tight')

plt.show()

