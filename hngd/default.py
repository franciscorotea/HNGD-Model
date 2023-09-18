# Default model parameters for the HNGD model

# Ref.:

# Passelaigue F, Lacroix E, Pastore G and Motta A T, "Implementation and 
# validation of the Hydride Nucleation-Growth-Dissolution (HNGD) model in
# BISON", Journal of Nuclear Materials 544 (2021), pp. 152683. 

# https://doi.org/10.1016/j.jnucmat.2020.152683

from numba.core import types
from numba.typed import Dict

# Use default model parameters, store them in a numba-compatible dictionary.

model_parameters = Dict.empty(key_type=types.unicode_type, 
                              value_type=types.float64)

model_parameters['TSSd0'] = 102000 # Pre-exponential term for TSSd [wt pmm]
model_parameters['Q_TSSd'] = 35458.7 # Activation energy for TSSd [J mol-1]
model_parameters['TSSp0'] = 30853.3 # Pre-exponential term for TSSp [wt ppm]
model_parameters['Q_TSSp'] = 25249.6 # Activation energy for TSSp [J mol-1]
model_parameters['D0'] = 1.08e-6 # Pre-exponential term of the diffusion coefficient [m2 s-1]
model_parameters['Ed'] = 0.46 # Activation energy of the diffusion coefficient [eV at-1]
model_parameters['Q_star'] = 25500 # Heat of transport of hydrogen in zirconium [J mol-1]
model_parameters['V_star'] = 1.7e-6 # Partial molar volume of hydrogen in the zirconium matrix [m3 mol-1]. CHECK > NORTHWOOD AND KOSASIH USE 6e-7 APPROX
model_parameters['Eg'] = 0.9 # Activation energy for diffusion driven growth [eV at-1]
model_parameters['avrami_parameter'] = 2.5 # Avrami parameter for platelets [/]

# Coefficients of the 3rd degree polynomial used to express the formation 
# energy of delta-hydrides

model_parameters['Eth0'] = 0.5655 # [eV at-1]
model_parameters['Eth1'] = 4.0e-4 # [eV at-1 K-1]
model_parameters['Eth2'] = 2.0e-7 # [eV at-1 K-2]
model_parameters['Eth3'] = 3.0e-10 # [eV at-1 K-3]

# Pre-exponential factors for growth, nucleation and dissolution kinetics

model_parameters['K_G_mob0'] = 5.53e5 # Pre-exponential factor for diffusion-controlled growth [s-1]
model_parameters['K_G_th0'] = 1.6e-5 # Pre-exponential factor for reaction-controlled growth [s-1]
model_parameters['K_N_0'] = 2.75e-5 # Pre-exponential factor for nucleation kinetics [s-1]
model_parameters['K_D_0'] = 1110.13 # Pre-exponential factor for dissolution kinetics [s-1]

# mHNGD model parameters

model_parameters['tau'] = 1e4 # s
model_parameters['g'] = 120 # wt.ppm
model_parameters['delta'] = 1.05 # /