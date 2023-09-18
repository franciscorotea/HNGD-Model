import matplotlib.pyplot as plt
import numpy as np
import hngd.default
import hngd.run

"""
Generate plots in Fig. 5.2 of Lacroix E., 'Modeling zirconium hydride 
precipitation and dissolution in zirconium alloys', PhD Thesis, The 
Pennsylvania State University p. 96 (2019).

This plot is used as a benchmark to verify precipitation and dissolution 
kinetics. Experimental data is taken from Lacroix E., Motta A. T, and Almer 
J. D., 'Experimental determination of zirconium hydride precipitation and 
dissolution in zirconium alloy', Journal of Nuclear Materials 509, pp. 162-167 
(2018).

"""

# Input data for simulation.

experiment_data = (0.025, 0, 0, 1, 0)

simulation_data = ('lacroix_experiment', 5, 2, 20, 0.001, 'mHNGD')

temperature_data = np.array([[np.nan,   0, experiment_data[0]],
                             [     0,  20,                 20],
                             [    33, 425,                425],
                             [    43, 425,                425],
                             [    65, 170,                170],
                             [    75, 170,                170],
                             [    90, 385,                385],
                             [   100, 385,                385],
                             [   110, 285,                285],
                             [   130, 285,                285],
                             [   140, 425,                425],
                             [   150, 425,                425],
                             [   170, 180,                180]])

initial_hydrogen_data = np.array([[  0, experiment_data[0]],
                                  [255,               255]])

model_parameters = hngd.default.model_parameters

#%% Run simulation.

hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, 
                    simulation_data, model_parameters)

#%% Plot

# Load experimental data and results from simulation.

exp_time, exp_Css = np.loadtxt('Data/lacroix_experimental_data.txt', skiprows=1, unpack=True)
sim_time, sim_temperature, _, _, sim_Css, _ = np.loadtxt(f'Outputs/{simulation_data[0]}_node_1_output_file.txt', 
                                                         skiprows=2, unpack=True)

# Plot.

fig, ax1 = plt.subplots(figsize=(8, 4))

fig.canvas.manager.set_window_title('Lacroix experiment')

ax1.scatter(exp_time, exp_Css, label='Experimental', facecolors='none', edgecolors='#1f77b4', s=15)
ax1.plot(sim_time, sim_Css, label='Simulated')

ax1.annotate('', [1500,250], [600,250], arrowprops={'arrowstyle':'<-', 'color':'#1f77b4'})

ax1.set_xlim([0, 450])
ax1.set_ylim([0, 300])

ax1.set_ylabel('Css [wt.ppm]', color='#1f77b4')
ax1.set_xlabel('Time [s]')

ax1.tick_params(axis='y', labelcolor='#1f77b4')

ax1.legend()

ax2 = ax1.twinx()

ax2.plot(sim_time, sim_temperature-273, color='#ff7f0e', linestyle='dashed')

ax1.annotate('', [9700,180], [10150,180], arrowprops={'arrowstyle':'<-', 'color':'#ff7f0e'})

ax2.set_ylabel('Temperature [°C]', color='#ff7f0e')
ax2.tick_params(axis='y', labelcolor='#ff7f0e')
ax2.set_xlim([sim_time[0], sim_time[-1]])
ax2.set_ylim([0, 550])

fig.tight_layout()

fig.savefig(f'Outputs/{simulation_data[0]}.png', bbox_inches='tight')

plt.show()
