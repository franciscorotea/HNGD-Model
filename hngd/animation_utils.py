"""A module with helper functions used to do matplolib animations."""

import numpy as np
import matplotlib.pyplot as plt

from .hngd import calc_terminal_solid_solubility

def init_animate(lines):
    """Initialize the animation with empty data."""
    for line in lines:
        line.set_data([], [])
    return lines

def animate(simulate, lines, texts, node, t_max, 
            time_evol, temp_evol, hydr_evol):
    """Main animation function. Values for each variable of interest are
    obtained at each output timestep, and each line is updated with these
    values to perform the animation."""
    
    # Unpack variables of interest to plot.
    
    Css, Cp, TSSp, TSSd, growth_fraction, nucleation_fraction, temperature, t, test_t, x = simulate
    
    # Unpack all plot lines and texts.
    
    Css_profile, Cp_profile, TSSp_profile, TSSd_profile, temperature_profile, node_temperature_profile, node_Css_profile, node_temperature_evolution, node_Css_evolution, temperature_evolution, Css_evolution = lines

    text_time, text_nucleation, text_growth = texts
    
    # Change of units.
    
    x = x*1e3 # m to mm
    temperature = temperature-273 # K to °C
    
    # Update data.
    
    Css_profile.set_data(x, Css)
    Cp_profile.set_data(x, Cp)
    TSSp_profile.set_data(x, TSSp)
    TSSd_profile.set_data(x, TSSd)
    temperature_profile.set_data(x, temperature)
    
    node_temperature_profile.set_data([x[node]], [temperature[node]])
    node_Css_profile.set_data([x[node]], [Css[node]])
    node_temperature_evolution.set_data([t], [temperature[node]])
    node_Css_evolution.set_data([temperature[node]], [Css[node]])
    
    text_time.set_text(f'time = {np.round(t,1)} / {t_max} s')
    text_nucleation.set_text(f'nucleation = {np.round(nucleation_fraction*100,0)}%')
    text_growth.set_text(f'growth = {np.round(growth_fraction*100,0)}%')
    
    time_evol.append(t)
    temp_evol.append(temperature[node])
    hydr_evol.append(Css[node])
    
    temperature_evolution.set_data(time_evol, temp_evol)
    Css_evolution.set_data(temp_evol, hydr_evol)

    return Css_profile, Cp_profile, TSSp_profile, TSSd_profile, temperature_profile, node_temperature_profile, node_Css_profile, node_temperature_evolution, node_Css_evolution, temperature_evolution, Css_evolution, text_time, text_nucleation, text_growth
    
def set_up_figure(temperature_data, initial_hydrogen_data, exp_data, sim_data, 
                  model_parameters, node, ylim):
    """Set up the figure and the axes used to do the animation."""
    
    # Parse input data.
    
    X = sim_data[1]
    Nx = exp_data[0]
    liner_width = exp_data[4]
    liner_solubility_factor = exp_data[3]
    t_max = temperature_data[1:, 0][-1]*60
    
    dx = Nx/X
    x = np.linspace(-dx, Nx+dx, X+3)
    
    max_temperature = np.max(temperature_data[1:,1:])
    min_temperature = np.min(temperature_data[1:,1:])
    
    max_hydr = np.max(initial_hydrogen_data[1])
    
    # Create figure and axes.
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18,9))
    
    fig.canvas.manager.set_window_title('Animation of temperature/hydrogen evolution')
    
    # Set limits.
    
    ax1.set_xlim([0, Nx*1e3])
    ax2.set_xlim([0, Nx*1e3])
    ax3.set_xlim([min_temperature, max_temperature*1.2])
    ax4.set_xlim([0, t_max])
    
    if ylim is not None:
        ax1.set_ylim([0, ylim])
    else:
        ax1.set_ylim([0, max_hydr*1.5])

    ax2.set_ylim([min_temperature, max_temperature*1.2])
    ax3.set_ylim([0, max_hydr*1.5])
    ax4.set_ylim([min_temperature, max_temperature*1.2])
    
    # Set labels.
    
    ax1.set_xlabel('Length [mm]')
    ax1.set_ylabel('Hydrogen concentration [wt.ppm]')

    ax2.set_xlabel('Length [mm]')
    ax2.set_ylabel('Temperature [°C]')
    
    ax3.set_xlabel('Temperature [°C]')
    ax3.set_ylabel('Hydrogen in solid solution concentration Css [wt.ppm]')
    
    ax4.set_xlabel('Time [s]')
    ax4.set_ylabel('Temperature [°C]')
    
    # Set titles.
    
    ax1.set_title('Hydrogen profile')
    ax2.set_title('Temperature profile')
    ax3.set_title(f'Evolution of hydrogen in solid solution at a single node N = {node} ({np.round(x[node]*1e3, 1)} mm)')
    ax4.set_title(f'Evolution of temperature at a single node N = {node} ({np.round(x[node]*1e3, 1)} mm)')
    
    # Plot TSSd and TSSp vs temperature (to calculate TSS, temperature needs to
    # be in K).
    
    TSS_temperature = np.linspace(1, max_temperature*1.2+273, 100)
    
    TSSd0 = model_parameters['TSSd0']
    Q_TSSd = model_parameters['Q_TSSd']
    TSSp0 = model_parameters['TSSp0']
    Q_TSSp = model_parameters['Q_TSSp']

    TSSp, TSSd = calc_terminal_solid_solubility(TSSp0, Q_TSSp, 
                                                TSSd0, Q_TSSd, 
                                                TSS_temperature, x, liner_width,
                                                liner_solubility_factor)
    
    ax3.plot(TSS_temperature-273, TSSp, ':k', label='TSSp')
    ax3.plot(TSS_temperature-273, TSSd, '--k', label='TSSd')
    
    ax3.legend(loc='upper left')
    
    fig.tight_layout()
    
    return fig, ((ax1, ax2), (ax3, ax4))

def set_up_lines(ax1, ax2, ax3, ax4):
    """Define all plot lines to be animated. Data will change at each time step,
    and all plots will be updated correspondingly in the animate function."""
    
    # ax1: Hydrogen profile
    
    Css_profile, = ax1.plot([], [], lw=2, ls='-', marker='o', ms=3, label='Css')
    Cp_profile, = ax1.plot([], [], lw=2, ls='-', marker='o', ms=3, label='Cp')
    TSSp_profile, = ax1.plot([], [], lw=1, ls='--', marker='s', ms=2, label='TSSp')
    TSSd_profile, = ax1.plot([], [], lw=1, ls=':', marker='s', ms=2, label='TSSd')
    node_Css_profile, = ax1.plot([], [], marker='x', mew=2, color='#9467bd')
    
    text_time = ax1.text(0.03, 0.97, '', horizontalalignment='left',
                         verticalalignment='top', transform=ax1.transAxes)
    
    text_nucleation = ax1.text(0.03, 0.92, '', horizontalalignment='left',
                               verticalalignment='top', transform=ax1.transAxes)
    
    text_growth = ax1.text(0.03, 0.875, '', horizontalalignment='left',
                           verticalalignment='top', transform=ax1.transAxes)
    
    ax1.legend(loc='upper right', ncol=2)
    
    # ax2: Temperature profile
    
    temperature_profile, = ax2.plot([], [], lw=2, ls='-', marker='o', ms=3)
    node_temperature_profile, = ax2.plot([], [], marker='x', mew=2, color='#9467bd')
    
    # ax3: Evolution of Css at a single node
    
    node_Css_evolution, = ax3.plot([], [], 'x', color='#9467bd')
    Css_evolution, = ax3.plot([], [], linestyle=':', marker='x', color='#9467bd', alpha=0.5)
    
    # ax4: Evolution of temperature at a single node
    
    node_temperature_evolution, = ax4.plot([], [], 'x', color='#9467bd')
    temperature_evolution, = ax4.plot([], [], linestyle=':', marker='x', color='#9467bd', alpha=0.5)
    
    lines = (Css_profile, Cp_profile, TSSp_profile, TSSd_profile, temperature_profile, 
             node_temperature_profile, node_Css_profile, node_temperature_evolution, node_Css_evolution, 
             temperature_evolution, Css_evolution)
    
    texts = (text_time, text_nucleation, text_growth)
    
    return lines, texts

def animate_temperature(run_thermal_history, lines, texts, t_max):
    """Animation function for thermal history."""
    
    # Unpack variables of interest to plot.
    
    t, x, temperature_profile, input_temperature_positions, temperatures_at_input_pos = run_thermal_history
    
    # Unpack all plot lines and texts.
    
    interp_temperature, input_temperature = lines

    text_time, = texts
    
    # Update data.
    
    interp_temperature.set_data(x, temperature_profile)
    input_temperature.set_data(input_temperature_positions, temperatures_at_input_pos)
    
    text_time.set_text(f'time = {np.round(t,1)} / {t_max} min')

    return interp_temperature, input_temperature, text_time