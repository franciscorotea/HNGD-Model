import HNGD.run
import HNGD.default
import experimental_data
import numpy as np
import matplotlib.pyplot as plt

# TODO > En la tesis de Florian Passelaigue (2020 - Hydride 
# nucleation-growth-dissolution model: Implementation in BISON) muestra
# barras de error para estos experimentos (no encontre de donde salen?).

"""
Generate plots for the experiments made by Kammenzind, published in Merlino 
J. T., 'Experiments in hydrogen distribution in thermal gradients calculated
using BISON', M. Eng. Thesis, The Pennsylvania State University (2019).

This script runs experiments from the Experiment 2 Data Set. Experimental data
can be found in Appendix A, Tables A15-A18. In these experiments, specimens 
were plated with a hydride rim on the hotter end, instead of being 
homogeneously distributed with hydrogen. In this case, the annealing 
temperature gradient did not increase linearly but had a colder side 
temperature of either 260°C or 316°C, a peak temperature of 371°C about 2.54 cm 
from the colder side, and a hotter side temperature at the opposite end of the 
specimen with a temperature between the colder side and the peak temperature. 
The anneals were conducted for periods of time ranging from 95 to 209 days. The 
specimens were then removed from heat, sectioned, and analyzed for total 
hydrogen content. 

Zirconium alloys used include alpha-annealed Zircaloy-4, Zircaloy-1%Niobium and
Zircaloy-2.5%Niobium.

This script takes ~ 60 minutes to run all 5 experiments. 

"""

slab_length = 3.810*1e-2
hydride_rim = 0.100*1e-2

experiment_data = (slab_length, 0, 0, 1, 0)

model_parameters = HNGD.default.model_parameters

# TSSp and TSSd of Zry-4 by Zanellato et al.

model_parameters['TSSd0'] = 510800
model_parameters['Q_TSSd'] = 45610
model_parameters['TSSp0'] = 66440
model_parameters['Q_TSSp'] = 29630

# Diffusion coefficient for Zry-4 by Kammenzind et al.

#model_parameters['D0'] = 0.8e-7
#model_parameters['Ed'] = 0.34

# Diffusion coefficient for Zry-4 by Kearns.

model_parameters['D0'] = 7.73e-7
model_parameters['Ed'] = 0.47

#%% Simulate

for i, (experiment, data) in enumerate(experimental_data.asymmetric_cases.items()):
    
    print(f'\n- Running experiment {experiment}. Simulation {i+1} out of {len(experimental_data.asymmetric_cases)}.')
    
    simulation_data = (experiment, 80, 2, 50000, 0.05, 'mHNGD')
    temperature_data = data['temperature_data']
    initial_hydrogen_data = np.array([[ 0, slab_length-hydride_rim, slab_length-hydride_rim+1e-12,               slab_length],
                                      [10,                      10,      data['initial_hydrogen'], data['initial_hydrogen']]])
    model_parameters['Q_star'] = data['heat_of_transport']
    
    HNGD.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, 
                        simulation_data, model_parameters, save_results=False, 
                        save_last_concentration_profile=True)
    
#%% Plot

for experiment, data in experimental_data.asymmetric_cases.items():
    
    # Load experimental and simulated data.
    
    pos_exp = data['midpoint_location']
    hyd_exp = data['hydrogen_concentration']
    
    pos_sim, temp_sim, Css_sim, Cp_sim = np.loadtxt(f'Outputs/{experiment}_last_concentration_profile.txt', 
                                                    skiprows=1, unpack=True)
    hyd_sim = Css_sim + Cp_sim
    
    # Plot.
    
    fig, ax1 = plt.subplots()
    
    fig.canvas.manager.set_window_title(f'Simulation of the Kammenzind asymmetric case {experiment}')
    
    ax1.scatter(pos_exp*1e2, hyd_exp, label='Experimental')
    ax1.plot(pos_sim*1e2, hyd_sim, label='Simulated')

    ax1.set_xlabel('Distance [cm]')
    ax1.set_ylabel('Hydrogen content [wt.ppm]', color='#1f77b4')
    
    ax1.set_xlim([0, slab_length*1e2])
    ax1.set_ylim(bottom=0)
    
    ax1.tick_params(axis='y', labelcolor='#1f77b4')
    
    ax1.legend(loc='upper left')
    
    ax2 = ax1.twinx()
    
    ax2.plot(pos_sim*1e2, temp_sim, color='#ff7f0e', linestyle='dashed')
    
    ax2.set_ylabel('Temperature [K]', color='#ff7f0e')
    
    ax2.tick_params(axis='y', labelcolor='#ff7f0e')
    
    fig.savefig(f'Outputs/asymmetric_{experiment}.png', bbox_inches='tight')
    
plt.show()
