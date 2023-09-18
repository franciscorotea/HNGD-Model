import hngd.run
import hngd.default
import hngd.hngd
import numpy as np
import matplotlib.pyplot as plt

"""

Generate plot for a TTT diagram such as the one shown in Fig. 3.12 of Lacroix 
E., 'Modeling zirconium hydride precipitation and dissolution in zirconium 
alloys', PhD Thesis, The Pennsylvania State University p. 64 (2019).

All hydrogen is present in solid solution at the start of the simulation (since 
the initial temperature is above the dissolution temperature). Then, the sample 
is cooled up to the desired soak temperature at some cooling rate. Once this 
temperature is achieved, it is mainteined long enough so that the completion 
target of the reaction is achieved (in general, 0.99). Schematically:

initial_temperature [°C]
        \
         \
          \
           \     slope = cooling_rate [°C/min]
            \
             \
              \
               \_____________________     soak
                                      temperature
      <------->                          [°C]
    cooling_time [min]

This script takes ~ 1 minute to run. 

"""

min_temperature = 250 # Minimum temperature to calculate the TTT diagram [°C] 
max_temperature = 350 # Maximum temperature to calculate the TTT diagram [°C] 
temperature_step = 5 # Temperature steps to calculate the TTT diagram [°C]
cooling_rate = 10 # Cooling rate used to achieve the desired temperature [°C min-1]
completion_target = 0.99 # Desired reaction completion state [/]

experiment_data = (0.025, 0, 0, 1, 0)
simulation_data = ('TTT experiment', 8, 2, 80, 0.01, 'HNGD')

C_0 = 205.0

hydrog_data = np.array([[  0, experiment_data[0]],
                        [C_0,               C_0]])

model_parameters = hngd.default.model_parameters

#%% Simulate TTT experiments

hngd.run.TTT(hydrog_data, experiment_data, simulation_data, model_parameters, 
             min_temperature, max_temperature, temperature_step, 
             cooling_rate, completion_target)

#%% Plot TTT diagram

fig, ax = plt.subplots(figsize=(8,4))

fig.canvas.manager.set_window_title('TTT Diagram') 

# Load simulation data from the file.

TTT_sim_time, TTT_sim_temp = np.loadtxt('Outputs/TTT experiment.txt', skiprows=1, unpack=True)

# Experimental data

TTT_exp_time = [13715.421, 11586.75, 10906.245, 8944.371, 11784.18] # s
TTT_exp_temp = [280, 290, 301, 305, 320] # °C

ax.semilogx(TTT_sim_time, TTT_sim_temp, ':x', label='Simulated')
ax.semilogx(TTT_exp_time, TTT_exp_temp, 'o', fillstyle='none', label='Experimental')

#ax.set_xlim([1e3, 1e5])

ax.legend()

ax.set_xlabel('Time [s]')
ax.set_ylabel('Temperature [°C]')

fig.tight_layout()

fig.savefig('Outputs/TTT_diagram.png', bbox_inches='tight')

#%% Plot Css vs temperature

fig, ax = plt.subplots(figsize=(8,4))

fig.canvas.manager.set_window_title('Evolution of Css') 

TSS_temperature = np.linspace(230, 450, 100) + 273

TSSd0 = model_parameters['TSSd0']
Q_TSSd = model_parameters['Q_TSSd']
TSSp0 = model_parameters['TSSp0']
Q_TSSp = model_parameters['Q_TSSp']

TSSp, TSSd = hngd.hngd.calc_terminal_solid_solubility(TSSp0, Q_TSSp, TSSd0, Q_TSSd, TSS_temperature)

ax.plot(TSS_temperature-273, TSSp, ':k', label='TSSP')
ax.plot(TSS_temperature-273, TSSd, '--k', label='TSSD')

for soak_temperature in TTT_sim_temp.astype(int):
    time, temperature, _, _, Css, _ = np.loadtxt(f'Outputs/TTT_{soak_temperature}C_node_1_output_file.txt', 
                                                 skiprows=2, unpack=True)
    ax.plot(temperature-273, Css, label=f'{soak_temperature}°C')

ax.set_xlim([230, 450])
ax.set_ylim([0, 330])

# Two labels: one for TSSp and TSSd (legend1), and one for all experiments (legend2).

lines = plt.gca().get_lines()

legend1 = ax.legend([lines[i] for i in range(0,2)], [lines[i].get_label() for i in range(0,2)], loc='lower right')
legend2 = ax.legend([lines[i] for i in range(2,len(TTT_sim_temp)+2)], [lines[i].get_label() for i in range(2,len(TTT_sim_temp)+2)], loc='upper center', ncol=6)

plt.gca().add_artist(legend1)

ax.set_xlabel('Temperature [°C]')
ax.set_ylabel('Hydrogen in solid solution Css [wt ppm]')

fig.tight_layout()

fig.savefig('Outputs/Css_vs_temperature.png', bbox_inches='tight')

plt.show()