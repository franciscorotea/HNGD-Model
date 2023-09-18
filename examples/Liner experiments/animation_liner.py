import numpy as np
import hngd.default
import hngd.run

solubility_factor = 0.3 # from 0.3 to 1.0

# Input data

liner_width = 150e-6 # [m]
slab_length = liner_width + 550e-6 # [m]
hydride_rim = slab_length - 30e-6 # [m]
X = 100 # number of nodes in the sample
C_0 = 375 * X / 5 # total hydrogen concentration [wt ppm]

initial_temperature = 20 # [°C]
cold_temperature = 345 # [°C]
hot_temperature = 360 # [°C]
heating_rate = 20 # [°C/min]
cooling_rate = 0.5 # [°C/min]

heating_time = np.abs(initial_temperature-hot_temperature)/heating_rate
holding_time = heating_time + 6*60 # [min]
cooling_time = holding_time + np.abs(hot_temperature-initial_temperature)/cooling_rate

temperature_data = np.array([[      np.nan,                   0,         slab_length],
                             [           0, initial_temperature, initial_temperature],
                             [heating_time,     hot_temperature,    cold_temperature],
                             [holding_time,     hot_temperature,    cold_temperature],
                             [cooling_time, initial_temperature, initial_temperature]])

initial_hydrogen_data = np.array([[0, hydride_rim, hydride_rim+1e-9, slab_length],
                                  [0,           0,              C_0,        C_0]])

experiment_data = (slab_length, 0, 0, solubility_factor, liner_width)
simulation_data = (f'{solubility_factor}', 100, 2, 300, 0.05, 'HNGD')

model_parameters = hngd.default.model_parameters

# Animate

anim = hngd.run.animation(temperature_data, initial_hydrogen_data, 
                          experiment_data, simulation_data,
                          model_parameters, interval=2)