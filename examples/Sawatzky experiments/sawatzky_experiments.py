import matplotlib.pyplot as plt
import numpy as np
import hngd.default
import hngd.run

# TODO > En la tesis de Florian Passelaigue (2020 - Hydride 
# nucleation-growth-dissolution model: Implementation in BISON) muestra
# barras de error para estos experimentos (no encontre de donde salen?).

#%% Experimental data

"""
Generate plots in Fig. 5.3 of Lacroix E., 'Modeling zirconium hydride 
precipitation and dissolution in zirconium alloys', PhD Thesis, The 
Pennsylvania State University p. 99 (2019).

This plots reproduce the experimental results obtained in Sawatzky A., 
'Hydrogen in Zircaloy-2: Its distribution and heat of transport', Journal of 
Nuclear Materials 2, pp. 321-328 (1960). Experiments consist on two samples
charged with different amounts of hydrogen, that were then annealed under
different temperature gradients for a given time. Then, samples were cut into
slices to analyze the hydrogen content as a function of the position.

"""

experiment_data = (0.025, 0, 0, 1, 0)

model_parameters = hngd.default.model_parameters

a = {'name': 'sawatzky_experiment_1',
     'C_0': 130,
     'thermal_treatment': np.array([[   np.nan,   0, experiment_data[0]],
                                    [        0, 130,                477],
                                    [ 34*60*24, 130,                477]]),
     'title': '$C_0$ = 130 wt.ppm\n34-day anneal\n138.8 K/cm temperature gradient\n(temperature difference: 130-477°C)', 
     'color': '#1f77b4',
     'exp_pos': [0.04733, 0.2295, 0.39917, 0.58298, 0.77618, 0.95021, 1.12361, 1.31138, 1.48774, 1.67647, 1.85096, 2.02019, 2.21227],
     'exp_H_content': [202.38, 181.7, 317.74, 553.96, 481.50, 66.18, 36.73, 28.87, 23.99, 19.33, 18.84, 13.85, 14.3]}

b = {'name': 'sawatzky_experiment_2',
     'C_0': 64,
     'thermal_treatment': np.array([[   np.nan,   0, experiment_data[0]],
                                    [        0, 157,                454],
                                    [ 41*60*24, 157,                454]]),
     'title': '$C_0$ = 64 wt.ppm\n41-day anneal\n118.8 K/cm temperature gradient\n(temperature difference: 157-454°C)', 
     'color': '#1f77b4',
     'exp_pos': [0.04082, 0.22292, 0.40063, 0.60131, 0.77964, 0.95672, 1.13268, 1.30014, 1.47197, 1.6732, 1.84441, 2.11913],
     'exp_H_content': [188.415, 207.386, 305.787, 111.264, 28.839, 23.732, 21.707, 17.509, 18.192, 14.267, 14.575, 9.011]}

experiments = (a, b)

#%% Simulation

for i, experiment in enumerate(experiments):
    
    print(f'\n- Running Sawatzky experiment {experiment["name"][-1]}. Simulation {i+1} out of {len(experiments)}.')
    
    simulation_data = (experiment['name'], 100, 2, 80, 0.01, 'mHNGD')
    temperature_data = experiment['thermal_treatment']
    initial_hydrogen_data = np.array([[                0, experiment_data[0]],
                                      [experiment['C_0'], experiment['C_0']]])
    
    hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, 
                        simulation_data, model_parameters, save_results=False,
                        save_last_concentration_profile=True)
    
#%% Plot

for experiment in experiments:
    
    position, temperature, Css, Cp = np.loadtxt(f"Outputs/{experiment['name']}_last_concentration_profile.txt", 
                                                skiprows=1, unpack=True)
    
    fig, ax1 = plt.subplots()
    
    fig.canvas.manager.set_window_title(f"Sawatzky experiment {experiment['name'][-1]}")
    
    ax1.plot(position*1e2, Css+Cp, color=experiment['color'], label='Simulated')
    
    ax1.scatter(experiment['exp_pos'], experiment['exp_H_content'], color=experiment['color'], label='Experimental')
    
    ax1.axhline(y=experiment['C_0'], linestyle='dotted', color='k')
    
    ax1.annotate(f'$C_0$ = {experiment["C_0"]} wt.ppm', xy=(1.5,experiment["C_0"]), 
                 xytext=(0,3), textcoords='offset points', va='bottom')
    
    ax1.set_title(experiment['title'])
    
    ax1.set_xlim([0, experiment_data[0]*1e2])
    
    ax1.set_xlabel('Distance from cold side [cm]')
    ax1.set_ylabel('Hydrogen content [wt.ppm]', color=experiment['color'])
    ax1.tick_params(axis='y', labelcolor=experiment['color'])
    
    ax1.legend(loc='upper center')
    
    ax2 = ax1.twinx()
    
    ax2.plot(position*1e2, temperature-273, color='#ff7f0e', linestyle='dashed')
    
    ax2.set_ylabel('Temperature [°C]', color='#ff7f0e')
    ax2.tick_params(axis='y', labelcolor='#ff7f0e')
    
    fig.tight_layout()
    
    fig.savefig(f'Outputs/{experiment["name"]}.png', bbox_inches='tight')
    
plt.show()