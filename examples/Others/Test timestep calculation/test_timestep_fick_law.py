import hngd.run
import hngd.default
import numpy as np
import math
import matplotlib.pyplot as plt

"""

Verification of the diffusion process: Fick's law. Compare the analytical
solution of a diffusion couple of a 1 m long Zircaloy-4 bar with 100 wt.ppm
for half its length and 0 wt.ppm for the other half. Temperature is kept
constant at 700 K. Results show hydrogen migration towards the right side,
along the concentration gradient.

This script takes ~5 minutes to run.

Ref.:
    
Lee C. and Lee, Y., 'Simulation of hydrogen diffusion along the axial direction 
in zirconium cladding tube during dry storage', Journal of Nuclear Materials 
579, pp. 154352 (2023). 

https://doi.org/10.1016/j.jnucmat.2023.154352

"""

def analytical_solution_diffusion_couple(x, t):
    erf_array = np.array([math.erf((x_i-0.5)/(2*np.sqrt(D*t))) for x_i in x])
    return (C_left + C_right)/2 - (C_left - C_right)/2 * erf_array

temperature = 427 # °C

k_B = 8.617333262e-5 # eV/K
D0 = 1.08e-6 # m2/s
Ed = 0.46 # eV

D = D0 * np.exp(-Ed/(k_B*(temperature+273))) # m2/s

C_left = 100 # wt.ppm
C_right = 0 # wt.ppm

time = np.array([10, 20, 30])*60 # minutes

x = np.linspace(0, 1, 1000) # m

experiment_data = (1, 0, 0, 1, 0)

initial_hydrogen_data = np.array([[     0,    0.5, 0.5+1e-9,       1],
                                  [C_left, C_left,  C_right, C_right]])

model_parameters = hngd.default.model_parameters

#%% Simulate

for t in time:
    
    simulation_data = (f'fick_{t}', 1000, 2, 80, 0.01, 'HNGD')
    temperature_data = np.array([[np.nan,           0,           1],
                                 [     0, temperature, temperature],
                                 [     t, temperature, temperature]])
    
    hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, simulation_data, 
                        model_parameters, save_results=False, save_last_concentration_profile=True)

#%% Plot

fig, ax = plt.subplots()

fig.canvas.manager.set_window_title("Verification of Fick's law")

ax.plot([0, 0.5, 0.5, 1.0], [100, 100, 0, 0], color='k', ls='--', label='Initial hydrogen profile')

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

for i, t in enumerate(time):
    
    # Analytic
    
    analytic_concentration = analytical_solution_diffusion_couple(x, t*60) # time from minutes to seconds
    
    ax.plot(x, analytic_concentration, color=colors[i], label=f'Analytical, after {int(t/60)} h')
    
    # Simulation
    
    x_sim, _, Css, Cp = np.loadtxt(f'Outputs/fick_{t}_last_concentration_profile.txt', 
                                   skiprows=1, unpack=True)
    y_sim = Css + Cp

    ax.scatter(x_sim[1::2], y_sim[1::2], edgecolors=colors[i], facecolors='white', label=f'Simulated, after {int(t/60)} h')
    
    ax.legend()

ax.set_xlabel('Position [m]')
ax.set_ylabel('Hydrogen concentration [wt.ppm]')
    
ax.set_xlim([0.4750, 0.5250])
ax.set_ylim([-10, 110])
ax.grid()

fig.savefig("Outputs/fick_verification.png", bbox_inches='tight')

plt.show()