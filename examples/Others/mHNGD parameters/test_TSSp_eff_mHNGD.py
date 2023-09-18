import numpy as np
import matplotlib.pyplot as plt
import hngd.hngd
import hngd.default_model_parameters

"""

Reproduce Fig 3 from Passelaigue F., Simon, P. A. and Motta, A. T., 'Predicting 
the hydride rim by improving the solubility limits in the Hydride 
Nucleation-Growth-Dissolution (HNGD) model', Journal of Nuclear Materials 558, 
pp. 153363 (2022).
        
"""

model_parameters = hngd.default_model_parameters.model_parameters

tau = 1e4
t_0 = 0.3

temperature_stamps = np.array([20, 650, 300, 700, 700])
time_stamps = np.array([0, 0.05, 0.10, 0.3, 3])*1e4

time = np.linspace(time_stamps[0], time_stamps[-1], 1000)
temperature = np.interp(time, time_stamps, temperature_stamps)

TSSd0 = model_parameters['TSSd0']
Q_TSSd = model_parameters['Q_TSSd']
TSSp0 = model_parameters['TSSp0']
Q_TSSp = model_parameters['Q_TSSp']

temp = np.zeros(2)
TSSp = np.zeros(len(time))
TSSd = np.zeros(len(time))
TSSp_eff = np.zeros(len(time))

t_0 = 0

temperature_change_flag = True

for i, t in enumerate(time):
    temp[1] = temperature[i]
    
    TSSp[i], TSSd[i] = hngd.hngd.calc_terminal_solid_solubility(TSSp0, Q_TSSp, TSSd0, Q_TSSd, temp[1])
    TSSp_eff[i] = TSSp[i]
    
    diff_temp = temp[1] - temp[0]
    
    if diff_temp == 0 and temperature_change_flag:
        temperature_change_flag = False
        t_0 = t
    
    if diff_temp == 0:
        TSSp_eff[i] = TSSd[i] + (TSSp[i]-TSSd[i])*np.exp(-(t-t_0)/tau)
    else:
        temperature_change_flag = True
        t_0 = 0
        
    temp[0] = temp[1]

fig, ax1 = plt.subplots(figsize=(6,4))

color = 'tab:red'
ax1.plot(time, TSSd, color=color, label='TSSd', ls=':')
ax1.plot(time, TSSp, color=color, label='TSSp', ls='--')
ax1.plot(time, TSSp_eff, color=color, label='TSSp$_{eff}$')
ax1.set_ylim([0, 450])
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Hydrogen content [wt.ppm]', color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.legend()
ax1.grid()

ax2 = ax1.twinx()

color = 'tab:blue'
ax2.plot(time, temperature,color=color, label='Temperature')
ax2.set_ylabel('Temperature [K]', color=color)
ax2.tick_params(axis='y', labelcolor=color)

ax2.set_xlim([0, time_stamps[-1]])
ax2.set_ylim([0, 900])
ax2.legend()

fig.tight_layout()

plt.savefig('test_TSSp_eff_mHNGD.png', bbox_inches='tight')

plt.show()

