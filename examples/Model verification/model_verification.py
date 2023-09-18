import matplotlib.pyplot as plt
import numpy as np
import hngd.default
import hngd.run
import hngd.hngd

"""
Generate plots in Fig. 5.1 of Lacroix E., 'Modeling zirconium hydride 
precipitation and dissolution in zirconium alloys', PhD Thesis, The 
Pennsylvania State University p. 94 (2019).

These plots verify that the partial dierential equations are solved as they 
should and that the system behaves qualitatively as it is expected.

This script takes < 1 minute to run. 

"""

experiment_data = (0.025, 0, 0, 1, 0)

initial_hydrogen_data = np.array([[  0, experiment_data[0]],
                                  [200,               200]])

model_parameters = hngd.default.model_parameters

b = {'name': 'isothermal_precipitation',
     'thermal_treatment': np.array([[np.nan,   0, experiment_data[0]],
                                    [     0,  20,                 20],
                                    [   430, 450,                450],
                                    [   580, 300,                300],
                                    [   780, 300,                300],
                                    [  1000,  20,                20]]),
     'title': 'Isothermal precipitation', 
     'color': '#1f77b4',
     'arrows': [[[300,40], [350,88]], [[415,190], [448,190]], [[400,210], [340,210]], [[305,197], [288,170]], 
                [[284,120], [284,75]], [[280,64], [235,60]], [[160,43], [86,15]]]}

c = {'name': 'precipitation_with_hydrides', 
     'thermal_treatment': np.array([[np.nan,   0, experiment_data[0]],
                                    [     0,  20,                 20],
                                    [   360, 380,                380],
                                    [   720,  20,                20]]),
     'title': 'Precipitation with existing hydrides', 
     'color': '#ff7f0e',
     'arrows': [[[300,40], [350,88]], [[330,145], [290,115]], [[160,45], [100,20]]]}

d = {'name': 'precipitation_during_heating', 
     'thermal_treatment': np.array([[np.nan,   0, experiment_data[0]],
                                    [     0,  20,                 20],
                                    [   430, 450,                450],
                                    [   630, 250,                250],
                                    [   705, 300,                300],
                                    [   975,  20,                20]]),
     'title': 'Precipitation during heating', 
     'color': '#2ca02c',
     'arrows': [[[300,40], [350,88]], [[415,190], [448,190]], [[400,210], [340,210]], [[300,195], [245,120]], 
                [[267,96], [297,76]], [[270,65], [226,56]], [[150,40], [60,12]]]}

experiments = (b, c, d)

#%% Simulate

for i, experiment in enumerate(experiments):
    
    print(f'\n- Running {experiment["title"].lower()} experiment. Simulation {i+1} out of {len(experiments)}.')
    
    simulation_data = (experiment['name'], 5, 2, 8, 0.01, 'HNGD')
    temperature_data = experiment['thermal_treatment']
    
    hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, 
                        simulation_data, model_parameters)

#%% Plot

TSS_temperature = np.linspace(0, 500, 100) + 273

TSSd0 = model_parameters['TSSd0']
Q_TSSd = model_parameters['Q_TSSd']
TSSp0 = model_parameters['TSSp0']
Q_TSSp = model_parameters['Q_TSSp']

TSSp, TSSd = hngd.hngd.calc_terminal_solid_solubility(TSSp0, Q_TSSp, TSSd0, Q_TSSd, TSS_temperature)

for experiment in experiments:
    
    time, temperature, _, _, Css, _ = np.loadtxt(f"Outputs/{experiment['name']}_node_1_output_file.txt", 
                                                 skiprows=2, unpack=True)
    
    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    
    fig.canvas.manager.set_window_title(f'{experiment["title"]}')
    
    ax[0].plot(time, temperature-273, color=experiment['color'])
    
    ax[0].set_xlabel('Time [s]')
    ax[0].set_ylabel('Temperature [°C]')
    
    ax[0].set_xlim([time[0], time[-1]])
    ax[0].set_ylim([0, 500])
    
    ax[0].set_title('Heat treatment')
    
    ax[1].plot(TSS_temperature-273, TSSp, ':k', label='TSSP')
    ax[1].plot(TSS_temperature-273, TSSd, '--k', label='TSSD')
    
    ax[1].plot(temperature-273, Css, label='Css', color=experiment['color'])
    
    for arrow_loc in experiment['arrows']:
        ax[1].annotate('', *arrow_loc, arrowprops={'arrowstyle': '<-', 
                                                   'color': experiment['color'], 
                                                   'alpha': 0.5, 
                                                   'linewidth': 0.8})
    
    ax[1].set_xlabel('Temperature [°C]')
    ax[1].set_ylabel('Css [wt.ppm]')
    
    ax[1].set_xlim([0, 450])
    ax[1].set_ylim([0, 250])
    
    ax[1].set_title(experiment['title'])
    
    ax[1].legend()
    
    fig.tight_layout()
    
    fig.savefig(f"Outputs/{experiment['name']}.png", bbox_inches='tight')
    
plt.show()