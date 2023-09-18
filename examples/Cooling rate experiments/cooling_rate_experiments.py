import hngd.run
import hngd.default
import hngd.hngd
import numpy as np
import matplotlib.pyplot as plt

"""

Generate plot to simulate the influence of cooling rate on hydrogen behaviour, 
such as the ones shown in Fig. 5.4 and 5.5 of Lacroix E., 'Modeling zirconium 
hydride precipitation and dissolution in zirconium alloys', PhD Thesis, The 
Pennsylvania State University p. 102 (2019).

It is found that the fraction of hydrides generated from growth or nucleation
greatly depends on the cooling rate.

This script takes ~ 14 minutes to run (most of the time accounts for the 
0.01°C/min cooling rate experiment).

"""

# Cooling rates to test.

cooling_rates = np.array([0.01, 0.1, 0.2, 0.5, 0.8, 1.0, 2.0, 5.0, 8.0, 10.0, 
                          20.0, 50.0, 80.0, 100.0]) # [°C/min]

#%% Run simulations

# Input data

C_0 = 200 # wt.ppm
Nx = 0.025 # [m]

initial_temperature = 20 # [°C]
temperature = 450 # 450 for full dissolution or 380 for partial dissolution [°C]

model_parameters = hngd.default.model_parameters

initial_hydrogen_data = np.array([[  0,  Nx],
                                  [C_0, C_0]])

experiment_data = (Nx, 0, 0, 1, 0)

# Run simulations for different cooling rates.

for i, cooling_rate in enumerate(cooling_rates):
    
    print(f'\n- Cooling rate = {cooling_rate}°C/min. Simulation {i+1} out of {len(cooling_rates)}.')
    
    # Adjust the output criterion according to the length of the experiment.
    
    if 0.01 <= cooling_rate < 0.5:
        criterion = 0.05
        test_t_set = 180
        unit_time = 'd'        
    elif 0.5 <= cooling_rate < 2.0:
        criterion = 0.05
        test_t_set = 120
        unit_time = 'h'
    elif 2.0 <= cooling_rate <= 50.0:
        criterion = 0.01
        test_t_set = 10
        unit_time = 'm'
    else:
        criterion = 0.001
        test_t_set = 1
        unit_time = 's'
        
    simulation_data = (f'cooling_rate_{cooling_rate}', 8, 2, test_t_set, criterion, 'HNGD')
    
    heating_time = np.abs(initial_temperature-temperature)/cooling_rate
    holding_time = heating_time + 10
    cooling_time = holding_time + np.abs(temperature-initial_temperature)/cooling_rate

    temperature_data = np.array([[      np.nan,                   0,                  Nx],
                                 [           0, initial_temperature, initial_temperature],
                                 [heating_time,         temperature,         temperature],
                                 [holding_time,         temperature,         temperature],
                                 [cooling_time, initial_temperature, initial_temperature]])
    
    hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, 
                        simulation_data, model_parameters, unit_time=unit_time)
    #hngd.run.debugger(temperature_data, initial_hydrogen_data, experiment_data, 
    #                    simulation_data, model_parameters)
#%% Plot results

# Revisar > difiere un poco de la Figura 5.5 de la tesis para bajas tasas de 
# enfriamiento? Por ejemplo, para la tasa de 0.01 la fracción de precipitación
# debido al crecimiento deberia ser ~ 1 (por lo tanto, ~ 0 para la nucleacion)
# pero esta dando alrededor de ~ 0.8 (y ~ 0.2 para la nucleacion). Es raro
# porque la Figura 5.4 da practicamente igual para esa tasa de enfriamiento?

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

fig1.canvas.manager.set_window_title('Influence of cooling rate on hydrogen behaviour - Precipitation fraction vs cooling rate')
fig2.canvas.manager.set_window_title('Influence of cooling rate on hydrogen behaviour - Css vs temperature')

for cooling_rate in cooling_rates:
    
    time, temp, growth, nucleation, Css, _ = np.loadtxt(f'Outputs\cooling_rate_{cooling_rate}_node_1_output_file.txt', 
                                                            skiprows=2, unpack=True)
    
    if cooling_rate in [0.01, 0.1, 1.0, 10.0, 100.0]:
        ax2.plot(temp-273, Css, label=f'{cooling_rate}°C/min')
    
    fraction = np.column_stack((growth, nucleation))
    
    # Find the first index in which either nucleation or growth starts.
    
    idx = np.any(fraction != 0, axis=1).argmax(axis=0)
    
    # Filter the growth and nucleation fraction to take into account only 
    # values in which either nucleation or growth is present.
    
    fraction_filtered = fraction[idx:, :]
    
    # Find the dt for each fraction.
    
    dt_filtered = np.diff(time[idx-1:]-time[idx-1])
    
    # Weighted arithmetic mean.
    
    mean_fraction = np.sum(fraction_filtered * dt_filtered[:, np.newaxis], axis=0) / np.sum(dt_filtered)
    
    mean_growth = mean_fraction[0]
    mean_nucleation = mean_fraction[1]
    
    ax1.semilogx(cooling_rate, mean_nucleation, '-o', color='#1f77b4')
    ax1.semilogx(cooling_rate, mean_growth, '-x', color='#ff7f0e')

legend_1 = ax2.legend(loc='upper left')

# Plot TSSd and TSSp.

TSS_temperature = np.linspace(0, temperature, 100) + 273

TSSd0 = model_parameters['TSSd0']
Q_TSSd = model_parameters['Q_TSSd']
TSSp0 = model_parameters['TSSp0']
Q_TSSp = model_parameters['Q_TSSp']

TSSp, TSSd = hngd.hngd.calc_terminal_solid_solubility(TSSp0, Q_TSSp, TSSd0, Q_TSSd, TSS_temperature)

TSSp_plot, = ax2.plot(TSS_temperature-273, TSSp, ':k')
TSSd_plot, = ax2.plot(TSS_temperature-273, TSSd, '--k')

legend_2 = ax2.legend([TSSp_plot, TSSd_plot], ['TSSp', 'TSSd'], loc='lower right')

plt.gca().add_artist(legend_1)
plt.gca().add_artist(legend_2)

# Set axes labels and limits.

ax1.set_xlabel('Cooling rate [°C/min]')
ax1.set_ylabel('Precipitation fraction [/]')
ax2.set_xlabel('Temperature [°C]')
ax2.set_ylabel('Hydrogen in solid solution Css [wt.ppm]')

ax1.set_xlim([cooling_rates[0], cooling_rates[-1]])
ax1.set_ylim([0, 1])

ax2.set_xlim([0, temperature])
ax2.set_ylim([0, C_0*1.2])

ax1.legend(['Nucleation', 'Growth'], loc='upper left')
ax2.legend()

fig1.tight_layout()
fig2.tight_layout()

fig1.savefig('Outputs/precipitation_fraction_vs_cooling_rate.png', bbox_inches='tight')
fig2.savefig('Outputs/Css_vs_temperature.png', bbox_inches='tight')

plt.show()