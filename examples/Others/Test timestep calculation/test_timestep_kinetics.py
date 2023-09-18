import hngd.run
import hngd.default
import numpy as np
import matplotlib.pyplot as plt

# TODO > Hay algo raro en caso de crecimiento > en el paper dice que empieza con
# Css_init = 255-145 = 110 ppm, pero en el grafico se ve que empieza por encima
# de 175 ppm? Ademas, la solucion analitica no tiene la misma forma que se ve
# en el grafico del paper (la simulada si). 

# En el caso de la nucleacion, los resultados dan mas cercanos a la solucion
# analitica que en el paper (relacionado al dt calculado para el proceso de 
# nucleacion?)

"""

Verification of the time step used in the precipitation kinetics processes.

This script takes ~3 minutes to run.

Ref.:
    
Lee C. and Lee, Y., 'Simulation of hydrogen diffusion along the axial direction 
in zirconium cladding tube during dry storage', Journal of Nuclear Materials 
579, pp. 154352 (2023). 

https://doi.org/10.1016/j.jnucmat.2023.154352

"""

ppm_to_atomic_fraction = lambda x : x / (m_H * (x/m_H + (1e6-x)/m_Zr))

TSSp_fun = lambda T : TSSp0 * np.exp(-Q_TSSp/(R*T))
TSSd_fun = lambda T : TSSd0 * np.exp(-Q_TSSd/(R*T))

def calc_K_G(temperature, Cp, Css, TSSd):
    Eth = -Eth0 + Eth1*temperature - Eth2*temperature**2 + Eth3*temperature**3  
    TSSd_at = ppm_to_atomic_fraction(TSSd)
    Co_at = ppm_to_atomic_fraction(Cp + Css)
    Cp_at = ppm_to_atomic_fraction(Cp)
    C_delta = -9.93e-11*temperature**3 + 8.48e-8*temperature**2 - 5.73e-5*temperature + 0.623
    lvl_rule_Co = (Co_at - TSSd_at) / (C_delta - TSSd_at)
    lvl_rule_Cp = Cp_at / (C_delta - TSSd_at)
    Ko = K_G_mob0 * lvl_rule_Co * (1-lvl_rule_Cp)
    K1 = K_G_th0 * lvl_rule_Co * (1-lvl_rule_Cp) 
    K_mob = Ko * np.exp(-Eg/(k_B*temperature))
    K_ther = K1 * np.exp(-Eth/(k_B*temperature))
    K_G = 1 / (1/K_mob + 1/K_ther)
    return K_G    

def calc_K_D(temperature):
    K_D = K_D_0 * np.exp(-Ed/(k_B*temperature))
    return K_D

def calc_K_N(temperature, TSSd, Cp):
    Eth = -Eth0 + Eth1*temperature - Eth2*temperature**2 + Eth3*temperature**3
    TSSd_at = ppm_to_atomic_fraction(TSSd)
    Cp_at = ppm_to_atomic_fraction(Cp)
    C_delta = -9.93e-11*temperature**3 + 8.48e-8*temperature**2 - 5.73e-5*temperature + 0.623
    lvl_rule_Cp = Cp_at / (C_delta - TSSd_at)
    Kn0 = K_N_0 * (1 - lvl_rule_Cp)
    K_N = np.min([1, (17000-Cp)/17000*Kn0*np.exp(-Eth/(k_B*temperature))])
    return K_N

# Constants

R = 8.31446261815324
k_B = 8.617333262e-5
m_H = 1.008
m_Zr = 91.224

TSSd0 = 510800
Q_TSSd = 45610
TSSp0 = 66440
Q_TSSp = 29630
D0 = 1.08e-6
Ed = 0.46
Q_star = 25500
Eg = 0.9
p = 2.5

Eth0 = 0.5655
Eth1 = 4.0e-4
Eth2 = 2.0e-7
Eth3 = 3.0e-10

K_G_mob0 = 5.53e5
K_G_th0 = 1.6e-5
K_N_0 = 2.75e-5
K_D_0 = 1110.13

dissolution = {'initial_Css': 0,
               'initial_Cp': 550,
               'temperature': 600,
               'time': np.linspace(0, 35, 1000)}

nucleation = {'initial_Css': 550,
              'initial_Cp': 0,
              'temperature': 600,
              'time': np.linspace(0, 350, 1000)}

growth = {'initial_Css': 110,
          'initial_Cp': 145,
          'temperature': 603,
          'time': np.linspace(0, 7500, 1000)}

# Simulation parameters

experiment_data = (0.1, 0, 0, 1, 0)

initial_hydrogen_data = np.array([[  0,   0.1],
                                  [550,   550]])

model_parameters = hngd.default.model_parameters

model_parameters['TSSd0'] = 510800
model_parameters['Q_TSSd'] = 45610
model_parameters['TSSp0'] = 66440
model_parameters['Q_TSSp'] = 29630

#%% Create plot

fig, ax = plt.subplots(1, 3, figsize=(15,4))

fig.canvas.manager.set_window_title("Verification of the precipitation kinetics")

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

#%% Dissolution

# Simulation

initial_hydrogen_data = np.array([dissolution['initial_Css'], dissolution['initial_Cp']])
temperature_data = np.array([[np.nan,       0,     0.1],
                             [     0, 600-273, 600-273],
                             [ 35/60, 600-273, 600-273]])

for dt_div in [1, 2, 5, 10]:
    simulation_data = (f'dissolution_{dt_div}', 10, dt_div, 0, 0, 'HNGD')
    hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, simulation_data, 
                        model_parameters)

time1, _, _, _, Css1, _ = np.loadtxt('Outputs/dissolution_1_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)
time2, _, _, _, Css2, _ = np.loadtxt('Outputs/dissolution_2_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)
time5, _, _, _, Css5, _ = np.loadtxt('Outputs/dissolution_5_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)
time10, _, _, _, Css10, _ = np.loadtxt('Outputs/dissolution_10_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)

# Analytic

t = dissolution['time']
temperature = dissolution['temperature']
K_D = calc_K_D(temperature)
TSSd = TSSd_fun(temperature)
analytical_solution = TSSd * (1-np.exp(-K_D*t))

# Plot

ax[0].plot(time1, Css1, 'o:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=1')
ax[0].plot(time2, Css2, '^:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=2')
ax[0].plot(time5, Css5, 's:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=5')
ax[0].plot(time10, Css10, 'h:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=10')

ax[0].plot(t, analytical_solution, linewidth=2, color='k', label='Analytic')

ax[0].set_title('Dissolution')
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('Hydrogen in solid solution [wt.ppm]')
ax[0].legend()
ax[0].grid()

#%% Nucleation

# Simulation

initial_hydrogen_data = np.array([nucleation['initial_Css'], nucleation['initial_Cp']])
temperature_data = np.array([[np.nan,       0,     0.1],
                             [     0,  600-273, 600-273],
                             [ 350/60, 600-273, 600-273]])

for dt_div in [1, 2, 5, 10]:
    simulation_data = (f'nucleation_{dt_div}', 10, dt_div, 0, 0, 'HNGD')
    hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, simulation_data, 
                        model_parameters, activate_growth=False)

time1, _, _, _, Css1, _ = np.loadtxt('Outputs/nucleation_1_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)
time2, _, _, _, Css2, _ = np.loadtxt('Outputs/nucleation_2_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)
time5, _, _, _, Css5, _ = np.loadtxt('Outputs/nucleation_5_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)
time10, _, _, _, Css10, _ = np.loadtxt('Outputs/nucleation_10_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)

# Analytic

t = nucleation['time']
temperature = nucleation['temperature']
Cp = nucleation['initial_Cp']
Css = nucleation['initial_Css']
TSSp = TSSp_fun(temperature)
C_tot = Cp + Css
K_N = calc_K_N(temperature, TSSd, Cp)
analytical_solution = TSSp + (C_tot-TSSp) * np.exp(-K_N*t)

# Plot

ax[1].plot(time1, Css1, 'o:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=1')
ax[1].plot(time2, Css2, '^:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=2')
ax[1].plot(time5, Css5, 's:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=5')
ax[1].plot(time10, Css10, 'h:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=10')

ax[1].plot(t, analytical_solution, linewidth=2, color='k', label='Analytic')

ax[1].set_title('Nucleation')
ax[1].set_xlabel('Time [s]')
ax[1].set_ylabel('Hydrogen in solid solution [wt.ppm]')
ax[1].legend()
ax[1].grid()

#%% Growth

# Simulation

initial_hydrogen_data = np.array([growth['initial_Css'], growth['initial_Cp']])
temperature_data = np.array([[ np.nan,       0,      0.1],
                             [      0, 603-273, 603-273],
                             [7500/60, 603-273, 603-273]])

for dt_div in [1, 2, 5, 10]:
    simulation_data = (f'growth_{dt_div}', 10, dt_div, 0, 0, 'HNGD')
    hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, simulation_data, 
                        model_parameters, activate_nucleation=False)

time1, _, _, _, Css1, _ = np.loadtxt('Outputs/growth_1_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)
time2, _, _, _, Css2, _ = np.loadtxt('Outputs/growth_2_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)
time5, _, _, _, Css5, _ = np.loadtxt('Outputs/growth_5_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)
time10, _, _, _, Css10, _ = np.loadtxt('Outputs/growth_10_node_1_output_file.txt', 
                                     skiprows=2, unpack=True)

# Analytic

t = growth['time']
temperature = growth['temperature']
Cp = growth['initial_Cp']
Css = growth['initial_Css']
TSSd = TSSd_fun(temperature)
C_tot = Cp + Css
K_G = calc_K_G(temperature, Cp, Css, TSSd)
analytical_solution = TSSd + (C_tot-TSSd) * np.exp(-(K_G*t)**p)

# Plot

ax[2].plot(time1, Css1, 'o:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=1')
ax[2].plot(time2, Css2, '^:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=2')
ax[2].plot(time5, Css5, 's:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=5')
ax[2].plot(time10, Css10, 'h:', markerfacecolor='white', markersize=4, label='Numerical, dt_div=10')

ax[2].plot(t, analytical_solution, linewidth=2, color='k', label='Analytic')

ax[2].set_title('Growth')
ax[2].set_xlabel('Time [s]')
ax[2].set_ylabel('Hydrogen in solid solution [wt.ppm]')
ax[2].legend()
ax[2].grid()

fig.tight_layout()

fig.savefig("Outputs/kinetics_verification.png", bbox_inches='tight')

plt.show()