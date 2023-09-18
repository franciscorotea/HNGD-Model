import hngd.run
import hngd.default
import numpy as np
import matplotlib.pyplot as plt

# TODO > En el paper no dan la solucion analitica? Trate de hacer la que esta
# en el Appendix B de la tesis de Florian Passelaigue (2020 - Hydride 
# nucleation-growth-dissolution model: Implementation in BISON) pero no da.

# Hice las simulaciones con 6, 12 y 18 dias en lugar de horas como dice el
# paper (creo que se confundieron, sino da mucho menos difusion, usando dias
# en lugar de horas da casi igual a los resultados que presentan).

"""

Verification of the diffusion process: Soret effect. Compare the analytical
solution of a 10 cm long Zircaloy-4 bar with a uniform initial hydrogen
concentration of 50 wt.ppm subjected to a linear temperature profile of 
623K-673K. Results show hydrogen migration from the high temperature region
to the low temperature region.

This script takes ~5 minutes to run.

Ref.:
    
Lee C. and Lee, Y., 'Simulation of hydrogen diffusion along the axial direction 
in zirconium cladding tube during dry storage', Journal of Nuclear Materials 
579, pp. 154352 (2023). 

https://doi.org/10.1016/j.jnucmat.2023.154352

"""

def analytical_solution_soret_effect(x, t):
    
    a = -Q_star * L * C_0 / R
    b = np.sum(((temperature-temperature_left)/((temperature+273)*(temperature_left+273))) * dx)
    
    C_eq_0 = a * (1/b)
    
    C = C_eq_0 * np.exp(Q_star/R * (1/temperature - 1/temperature_left))
    
    return C

temperature_left = 350 # °C
temperature_right = 400 # °C

temperature = np.linspace(temperature_left, temperature_right, 100)

L = 0.1 # m
dx = L/100
C_0 = 50 # wt.ppm

Q_star = 25500 # J/mol
R = 8.31446261815324 # J mol-1 K-1
time = np.array([6, 12, 18])*24*60 # minutes

x = np.linspace(0, L, 100) # m

experiment_data = (L, 0, 0, 1, 0)

initial_hydrogen_data = np.array([[  0,  0.1],
                                  [C_0, C_0]])

model_parameters = hngd.default.model_parameters

model_parameters['TSSd0'] = 510800
model_parameters['Q_TSSd'] = 45610
model_parameters['TSSp0'] = 66440
model_parameters['Q_TSSp'] = 29630

#%% Simulate

for t in time:
    
    simulation_data = (f'soret_{t}', 10, 2, 80, 0.01, 'HNGD')
    temperature_data = np.array([[np.nan,                0,                 L],
                                 [     0, temperature_left, temperature_right],
                                 [     t, temperature_left, temperature_right]])
    
    hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, simulation_data, 
                        model_parameters, save_results=False, save_last_concentration_profile=True)

#%% Plot

fig, ax = plt.subplots()

fig.canvas.manager.set_window_title("Verification of Soret effect")

ax.plot([0, L], [C_0, C_0], color='k', ls='--', label='Initial hydrogen profile')

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

for i, t in enumerate(time):
    
    # Analytic
    
    analytic_concentration = analytical_solution_soret_effect(x, 0)
    
    ax.plot(x, analytic_concentration, color=colors[i], label=f'Analytical, after {int(t/(60*24))} days')
    
    # Simulation
    
    x_sim, _, Css, Cp = np.loadtxt(f'Outputs/soret_{t}_last_concentration_profile.txt', 
                                skiprows=1, unpack=True)
    y_sim = Css + Cp

    ax.plot(x_sim, y_sim, 'o:', markeredgecolor=colors[i], markerfacecolor='white', label=f'Simulated, after {int(t/(60*24))} days')
    
    ax.legend()

ax.set_xlabel('Position [m]')
ax.set_ylabel('Hydrogen concentration [wt.ppm]')
    
ax.set_xlim([0, L])
ax.set_ylim([45, 55])
ax.grid()

fig.savefig("Outputs/soret_verification.png", bbox_inches='tight')

plt.show()