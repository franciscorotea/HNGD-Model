import hngd.run
import hngd.default
import numpy as np

"""
A sample script to debug the hngd.hngd.simulate generator. Important: Set 
config.DISABLE_JIT = True in default_model_parameters to deactivate numba JIT
compilation.
"""
experiment_data = (0.1, 0, 0, 1, 0)
initial_hydrogen_data = np.array([550, 0])
temperature_data = np.array([[np.nan,       0,     0.1],
                             [     0, 600-273, 600-273],
                             [350/60, 600-273, 600-273]])

simulation_data = (f'dissolution_{1}', 10, 1, 0, 0, 'HNGD')

# =============================================================================
# experiment_data = (0.025, 0, 0, 1, 0)
# simulation_data = ('test', 70, 2, 80, 0.01, 'mhngd')
# 
# temperature_data = np.array([[np.nan,   0, experiment_data[0]],
#                              [     0,  20,                 20],
#                              [    38, 400,                400],
#                              [    48, 400,                400],
#                              [    86,  20,                20]])
# 
# initial_hydrogen_data = np.array([[  0, experiment_data[0]],
#                                   [205,               205]])
# =============================================================================

model_parameters = hngd.default.model_parameters

hngd.run.debugger(temperature_data, initial_hydrogen_data, experiment_data, 
                  simulation_data, model_parameters)
