import hngd.run
import hngd.default
import numpy as np
import matplotlib.pyplot as plt

"""

Generate plot to simulate the segregation of hydrides in liners, such as the 
one shown in Fig. 5.6 of Lacroix E., 'Modeling zirconium hydride precipitation 
and dissolution in zirconium alloys', PhD Thesis, The Pennsylvania State 
University p. 107 (2019).

A liner is a thin, protective layer that is applied to the inner surface of the
cladding. In the HNGD model, it is implemented through a local difference in 
solubility: liners are usually made of pure Zr and don't have much oxygen, 
causing the solubility in the liner to be lower than in the rest of the 
cladding. 

In this plot, different solubility factors are tested for typical reactor 
conditions. 

This script takes ~ 2 minutes to run. 

"""

# Input data

liner_width = 150e-6 # [m]
slab_length = liner_width + 550e-6 # [m]
hydride_rim = slab_length - 30e-6 # [m]
n_nodes = 100 # number of nodes in the sample

initial_temperature = 20 # [°C]
cold_temperature = 345 # [°C]
hot_temperature = 360 # [°C]
heating_rate = 20 # [°C/min]
cooling_rate = 0.5 # [°C/min]

heating_time = np.abs(initial_temperature-hot_temperature)/heating_rate # [min]
holding_time = heating_time + 6*60 # [min]
cooling_time = holding_time + np.abs(hot_temperature-initial_temperature)/cooling_rate # [min]

temperature_data = np.array([[      np.nan,                   0,         slab_length],
                             [           0, initial_temperature, initial_temperature],
                             [heating_time,    hot_temperature,     cold_temperature],
                             [holding_time,    hot_temperature,     cold_temperature],
                             [cooling_time, initial_temperature, initial_temperature]])

model_parameters = hngd.default.model_parameters

# All initial hydrogen is located in the last 30 um before the metal/oxide 
# interface, simulating the hydride rim tipically observed in nuclear reactors.

# Since the initial concentration of 375 wt ppm given in the thesis is not 
# homogeneously distributed across the sample, it has to be multiplied by a
# factor of (total number of nodes/number of nodes in the hydride rim).

C_0 = 375 * n_nodes / 5

initial_hydrogen_data = np.array([[0, hydride_rim, hydride_rim+1e-9, slab_length],
                                  [0,           0,              C_0,        C_0]])

# Run simulations for solubility factors from 0.3 to 1.0 in steps of 0.1.

solubility_factors = np.around(np.arange(0.3, 1.1, 0.1), 1)

for i, solubility_factor in enumerate(solubility_factors):
    
    print(f'\n- Solubility factor = {solubility_factor}. Simulation {i+1} out of {len(solubility_factors)}.')
    
    experiment_data = (slab_length, 0, 0, solubility_factor, liner_width)
    simulation_data = (f'solubility_factor_{solubility_factor}', n_nodes, 2, 300, 0.01, 'HNGD')

    hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, 
                        simulation_data, model_parameters, save_results=False, 
                        save_last_concentration_profile=True)
    
#%% Plot

# Si bien se observa un comportamiento similar al del grafico de la tesis de 
# Lacroix, no dan exactamente igual, en particular el contenido de hidrogeno
# total (alrededor de 10000 ppm en el rim en la tesis, y aca alrededor de 7500). 

# Disminuyendo un poco C_0 (de 7500 a 7000 wt ppm) se logra que todos 
# los hidruros pasen del rim al liner para sf = 0.3 (como se muestra en la 
# tesis).

fig, ax = plt.subplots(figsize=(8,4))

fig.canvas.manager.set_window_title('Hydrogen distribution for different solubility factors')

for solubility_factor in solubility_factors:

    position, temperature, Css, Cp = np.loadtxt(f'Outputs/solubility_factor_{solubility_factor}_last_concentration_profile.txt', 
                            skiprows=1, unpack=True)

    ax.plot(position*1e6, Css+Cp, label=f'sf = {solubility_factor}')

ax.axvline(x=liner_width*1e6, ls='--', color='k')

ax.text(0.10, 0.55, 'Liner', horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes)

ax.text(0.60, 0.55, 'Cladding', horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes)

ax.set_xlabel('Distance [μm]')
ax.set_ylabel('Hydrogen content [wt.ppm]')

ax.set_xlim([0, slab_length*1e6])
ax.set_ylim(bottom=0)

ax.legend(ncol=2, loc='upper center', title='solubility factor (sf)')

fig.tight_layout()

fig.savefig('Outputs/solubility_factors.png', bbox_inches='tight')

plt.show()