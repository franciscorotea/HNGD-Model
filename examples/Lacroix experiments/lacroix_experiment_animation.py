import numpy as np
import hngd.default
import hngd.run

# Input data for simulation.

experiment_data = (0.025, 0, 0, 1, 0)

simulation_data = ('lacroix_experiment', 5, 2, 20, 0.01, 'mHNGD')

temperature_data = np.array([[np.nan,   0, experiment_data[0]],
                             [     0,  20,                 20],
                             [    33, 425,                425],
                             [    43, 425,                425],
                             [    65, 170,                170],
                             [    75, 170,                170],
                             [    90, 385,                385],
                             [   100, 385,                385],
                             [   110, 285,                285],
                             [   130, 285,                285],
                             [   140, 425,                425],
                             [   150, 425,                425],
                             [   170, 180,                180]])

initial_hydrogen_data = np.array([[  0, experiment_data[0]],
                                  [255,               255]])

model_parameters = hngd.default.model_parameters

# Animate.

anim = hngd.run.animation(temperature_data, initial_hydrogen_data, 
                          experiment_data, simulation_data, model_parameters,
                          interval=50)