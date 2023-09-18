import hngd.run
import hngd.default
import numpy as np

experiment_data = (0.025, 0, 0, 1, 0)

simulation_data = ('test', 70, 2, 80, 0.01, 'mHNGD')

temperature_data = np.array([[np.nan,   0, experiment_data[0]],
                             [     0,  20,                 20],
                             [    38, 400,                400],
                             [    48, 400,                400],
                             [    86,  20,                20]])

initial_hydrogen_data = np.array([[  0, experiment_data[0]],
                                  [205,               205]])

model_parameters = hngd.default.model_parameters

hngd.run.simulation(temperature_data, initial_hydrogen_data, experiment_data, 
                    simulation_data, model_parameters)