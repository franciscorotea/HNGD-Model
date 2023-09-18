import numpy as np
import hngd.default
import hngd.run

# Input data.

cooling_rate = 10.0 # [°C/min]

C_0 = 200 # [wt ppm]
Nx = 0.025 # [m]

initial_temperature = 20 # [°C]
temperature = 450 # 450 for full dissolution or 380 for partial dissolution [°C]

model_parameters = hngd.default.model_parameters

initial_hydrogen_data = np.array([[  0,  Nx],
                                  [C_0, C_0]])

experiment_data = (Nx, 0, 0, 1, 0)
    
simulation_data = (f'example_cooling_rate_{cooling_rate}', 8, 2, 80, 0.01, 'HNGD')

heating_time = np.abs(initial_temperature-temperature)/cooling_rate
holding_time = heating_time + 10
cooling_time = holding_time + np.abs(temperature-initial_temperature)/cooling_rate

temperature_data = np.array([[      np.nan,                   0,                  Nx],
                             [           0, initial_temperature, initial_temperature],
                             [heating_time,         temperature,         temperature],
                             [holding_time,         temperature,         temperature],
                             [cooling_time, initial_temperature, initial_temperature]])

# Run animation.

anim = hngd.run.animation(temperature_data, initial_hydrogen_data, 
                          experiment_data, simulation_data,
                          model_parameters, interval=2)