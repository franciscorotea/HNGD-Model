import numpy as np
import hngd.default
import hngd.run

# Input data for simulation.

experiment = '1' # 1 or 2

experiment_data = (0.025, 0, 0, 1, 0)

data = {'1': {'C_0': 130,
              'ylim': 600,
              'thermal_treatment': np.array([[   np.nan,   0, experiment_data[0]],
                                             [        0, 130,                477],
                                             [ 34*60*24, 130,                477]])},
        '2': {'C_0': 130,
              'ylim': 600,
              'thermal_treatment': np.array([[   np.nan,   0, experiment_data[0]],
                                             [        0, 130,                477],
                                             [ 34*60*24, 130,                477]])}}

simulation_data = (experiment, 100, 2, 1000, 0.05, 'mHNGD')

initial_hydrogen_data = np.array([[                      0,      experiment_data[0]],
                                  [data[experiment]['C_0'], data[experiment]['C_0']]])

model_parameters = hngd.default.model_parameters

# Animate.

anim = hngd.run.animation(data[experiment]['thermal_treatment'], initial_hydrogen_data, 
                          experiment_data, simulation_data, model_parameters,
                          interval=2, ylim=data[experiment]['ylim'])