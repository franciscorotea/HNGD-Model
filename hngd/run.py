import os
import pdb

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation

from functools import partial
from tqdm import tqdm

from .hngd import simulate, run_thermal_history
from .io_utils import open_results_file, write_results, close_results_file
from .animation_utils import set_up_figure, set_up_lines, animate, init_animate, animate_temperature

progress_bar_format = '{desc}: {percentage:.1f}%|{bar}| {n:.1f}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]'

def simulation(temperature_data, initial_hydrogen_data, experiment_data, 
               simulation_data, model_parameters, node=1, 
               save_results=True, save_last_concentration_profile=False,
               unit_time='s', unit_temperature='K', activate_nucleation=True, 
               activate_growth=True, activate_dissolution=True):
    
    t_max = temperature_data[1:, 0][-1]*60

    print('\nStarting simulation...\n')
    
    # Set progress bar.
    
    pbar = tqdm(total=t_max, bar_format=progress_bar_format)
    
    # Prepare output file.
    
    if save_results:
        
        if not os.path.exists('Outputs'):
            os.makedirs('Outputs')
            
        output_file = open_results_file(f'Outputs/{simulation_data[0]}_node_{node}_output_file.txt', 
                                        node, unit_time, unit_temperature)
    
    # Run simulation.
    
    for sim_step in simulate(temperature_data, initial_hydrogen_data, experiment_data, simulation_data, model_parameters, 
                             activate_nucleation, activate_growth, activate_dissolution):
        
        Css, Cp, TSSp, TSSd, growth_fraction, nucleation_fraction, temperature_profile, t, dt, x = sim_step
        
        if save_results:
            write_results(output_file, t, temperature_profile[node], 
                          growth_fraction, nucleation_fraction, 
                          Css[node], Cp[node], unit_time, unit_temperature)
            
        pbar.update(dt)
    
    if save_last_concentration_profile:
        
        if not os.path.exists('Outputs'):
            os.makedirs('Outputs')
            
        if unit_temperature=='K':
            result = np.column_stack((x[1:-1], temperature_profile[1:-1], Css[1:-1], Cp[1:-1]))
            np.savetxt(f'Outputs/{simulation_data[0]}_last_concentration_profile.txt', result, 
                       header=f'{"Position [m]":^24} {"Temperature [K]":^24} {"Css [wt.ppm]":^24} {"Cp [wt.ppm]":^24}', 
                       comments='')
        elif unit_temperature=='C':
            result = np.column_stack((x[1:-1], temperature_profile[1:-1]-273, Css[1:-1], Cp[1:-1]))
            np.savetxt(f'Outputs/{simulation_data[0]}_last_concentration_profile.txt', result, 
                       header=f'{"Position [m]":^24} {"Temperature [°C]":^24} {"Css [wt.ppm]":^24} {"Cp [wt.ppm]":^24}', 
                       comments='')
        
        print(f'\nSimulation completed, last concentration profile saved at {os.getcwd()}\Outputs\{simulation_data[0]}_last_concentration_profile.txt.')

    # Close progress bar and output file.

    if save_results:
        print(f'\nSimulation completed, results for node {node} saved at {os.getcwd()}\Outputs\{simulation_data[0]}_node_{node}_output_file.txt.')
        close_results_file(output_file)
        
    pbar.close()
    
def TTT(hydrogen_data, experiment_data, simulation_data, model_parameters, 
        min_temperature, max_temperature, temperature_step, cooling_rate, 
        completion_target):
    
    # Calculating the precipitation and dissolution temperatures. 
    
    C_0 = hydrogen_data[1][0]
    R = 8.31446261815324 # Gas constant [J mol-1 K-1]
    
    TSSd0 = model_parameters['TSSd0']
    Q_TSSd = model_parameters['Q_TSSd']
    TSSp0 = model_parameters['TSSp0']
    Q_TSSp = model_parameters['Q_TSSp']

    precipitation_temperature = Q_TSSp/(R*np.log(TSSp0/C_0)) - 273 # °C
    dissolution_temperature = Q_TSSd/(R*np.log(TSSd0/C_0)) - 273 # °C

    # Set the initial temperature of the experiment so that all hydrogen is
    # in solid solution.

    initial_temperature = dissolution_temperature + 10

    # Temperatures are calculated in steps of `temperature_step`.

    soak_temperatures = np.arange(min_temperature, max_temperature, temperature_step)

    # The TTT diagram can't be calculated above the precipitation temperature.

    soak_temperatures = soak_temperatures[soak_temperatures < precipitation_temperature]

    # Run TTT experiments.

    print('\nStarting TTT simulation...\n')

    TTT_sim_time, TTT_sim_temp = [], []

    # Set up progress bar.

    pbar = tqdm(total=np.max(soak_temperatures)-np.min(soak_temperatures), 
                 bar_format=progress_bar_format, position=0, leave=True)

    for soak_temperature in soak_temperatures:
        
        # Create thermal treatment for this temperature.
        
        cooling_time = np.abs(initial_temperature-soak_temperature)/cooling_rate
        
        temperature_data = np.array([[            np.nan,                   0,  experiment_data[0]],
                                     [                 0, initial_temperature, initial_temperature],
                                     [      cooling_time,    soak_temperature,    soak_temperature],
                                     [cooling_time+36000,    soak_temperature,    soak_temperature]])
        
        temperature_flag = True
        
        if not os.path.exists('Outputs'):
            os.makedirs('Outputs')
            
        output_file = open_results_file(filename=f'Outputs\TTT_{soak_temperature}C_node_1_output_file.txt', node=1)
        
        # Simulate
        
        for sim_step in simulate(temperature_data, hydrogen_data, experiment_data, simulation_data, model_parameters):
            
            Css, Cp, TSSp, TSSd, growth_fraction, nucleation_fraction, temperature, t, dt, x = sim_step
            
            # Write the results (used to plot afterwards). Note that, as there
            # is no temperature/hydrogen gradient, the temperature/hydrogen at
            # any node would be the same (because of the flat profile).
            
            write_results(output_file, t, temperature[1], growth_fraction, nucleation_fraction, Css[1], Cp[1])
            
            # Calculate the advancement of the reaction.
            
            reaction_advancement = Cp[1] / (Cp[1] + Css[1] - TSSd[1])
            
            # Find the time needed to reach the soak temperature.
            
            if temperature[1]-273 == soak_temperature and temperature_flag:
                temperature_flag = False
                t_offset = t
                
            # Find the time needed to complete the reaction (the time
            # needed to reach the target temperature must be removed!). 
            
            if reaction_advancement >= completion_target and temperature[1]-273 == soak_temperature:
                TTT_sim_time.append(t-t_offset)
                TTT_sim_temp.append(soak_temperature)
                close_results_file(output_file)
                break
            
        pbar.update(temperature_step)
        
    pbar.close()

    # Save results.
    
    if not os.path.exists('Outputs'):
        os.makedirs('Outputs')
        
    result = np.column_stack((TTT_sim_time, TTT_sim_temp))
    
    np.savetxt(f'Outputs/{simulation_data[0]}.txt', result, 
               header=f'{"Time [s]":^24} {"Temperature [°C]":^24}', comments='')
    
    print(f'\nSimulation completed, TTT results saved at {os.getcwd()}\Outputs\{simulation_data[0]}.txt.')

def animation(temperature_data, initial_hydrogen_data, experiment_data, simulation_data, 
              model_parameters, node=None, save_animation=False, interval=5, repeat=False, ylim=None):
    
    if not node:
        node = simulation_data[1]//2 + 1
    
    t_max = temperature_data[1:, 0][-1]*60
    
    # Create figure and set up titles, labels, limits, etc.

    fig, ((ax1, ax2), (ax3, ax4)) = set_up_figure(temperature_data, 
                                                  initial_hydrogen_data, 
                                                  experiment_data, simulation_data, 
                                                  model_parameters, 
                                                  node, ylim)

    # Create all lines to be animated in each axis.

    lines, texts = set_up_lines(ax1, ax2, ax3, ax4)
    
    # Empty lists to animate temperature and hydrogen temporal evolution.
    
    time_evol, temp_evol, hydr_evol = [], [], []

    # Animate

    anim = matplotlib.animation.FuncAnimation(fig=fig, 
                                              func=animate,
                                              init_func=partial(init_animate, lines), 
                                              frames=partial(simulate, temperature_data, initial_hydrogen_data, experiment_data, simulation_data, model_parameters)(),
                                              fargs=(lines, texts, node, t_max, time_evol, temp_evol, hydr_evol),
                                              interval=interval,
                                              repeat=repeat, 
                                              blit=True,
                                              cache_frame_data=False)
    
    # TODO > The animation runs twice when save_animation=True, once for
    # saving and once for reproducing the animation. Do it at the same time?
    # Also, when saving animation xdata, ydata should be reset to empty lists?
    
    if save_animation:
        
        progress_callback = lambda i, n: print(f'Saving frame {i}')
        
        if not os.path.exists('Outputs'):
            os.makedirs('Outputs')
        
        anim.save(f'Outputs\{simulation_data[0]}_animation.gif', 
                  progress_callback=progress_callback)
        
        print(f'\nAnimation saved at {os.getcwd()}\Outputs\{simulation_data[0]}_animation.gif.')
        
    plt.show()
    
    return anim

def thermal_history_animation(temperature_data, n_nodes=10, dt=0.1, 
                              save_animation=False, interval=5, repeat=False):
    
    # Parse input data

    input_temperature_positions = temperature_data[0, 1:]
    input_time_stamps = temperature_data[1:, 0]
    input_temperature = temperature_data[1:, 1:]
    
    max_temperature = np.max(temperature_data[1:,1:])
    min_temperature = np.min(temperature_data[1:,1:])
    
    t_max = input_time_stamps[-1]
    
    # Create figure and axes.
    
    fig, ax = plt.subplots(figsize=(7,4))
    
    fig.canvas.manager.set_window_title('Animation of thermal history')
    
    ax.set_xlim([0, input_temperature_positions[-1]])
    ax.set_ylim([min_temperature, max_temperature*1.2])
    
    ax.set_xlabel('Length [m]')
    ax.set_ylabel('Temperature [°C]')

    # Create lines and text to be animated.
    
    interp_temperature, = ax.plot([], [], ls=':')
    input_temperature, = ax.plot([], [], marker='o', linestyle = 'None', mew=1, mfc='white', mec='#1f77b4')
    
    lines = (interp_temperature, input_temperature)
    
    text_time = ax.text(0.03, 0.97, '', horizontalalignment='left',
                        verticalalignment='top', transform=ax.transAxes)
    
    texts = (text_time,)
    
    # Animate

    anim = matplotlib.animation.FuncAnimation(fig=fig, 
                                              func=animate_temperature,
                                              init_func=partial(init_animate, lines), 
                                              frames=partial(run_thermal_history, temperature_data, n_nodes, dt)(),
                                              fargs=(lines, texts, t_max),
                                              interval=interval,
                                              repeat=repeat, 
                                              blit=True,
                                              cache_frame_data=False,
                                              save_count=int(t_max/dt))
    
    if save_animation:
        
        progress_callback = lambda i, n: print(f'Saving frame {i}')
        
        if not os.path.exists('Outputs'):
            os.makedirs('Outputs')
        
        anim.save('thermal_history_animation.gif',
                  progress_callback=progress_callback)
        
        print(f'\nAnimation of thermal history saved at {os.getcwd()}\Outputs\thermal_history_animation.gif.')
        
    fig.tight_layout()
    
    plt.show()
    
    return anim

def debugger(temperature_data, initial_hydrogen_data, experiment_data, 
             simulation_data, model_parameters):
    
    x = simulate(temperature_data, initial_hydrogen_data, experiment_data, 
                 simulation_data, model_parameters)
    
    results = []
    
    try:
        while True:
                result = pdb.runcall(next, x)
                results.append(result)
    except StopIteration:
        return results