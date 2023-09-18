import hngd.run
import hngd.default
import experimental_data
import numpy as np

# Choose an experiment from the data base. Valid options for linear cases are:
# A26a, A27, A45, A46, A26b, C13, A09a, A10a, A09b, A10b, A31, C14a, B14, 
# C14b, A11a, A12a, B15, C15a, A12b, C15b, A11b, C18, A13a, A14, C16a, B16, 
# A13, C16b.

experiment = 'A26a'

# Input parameters.

slab_length = 2.54*1e-2

data = experimental_data.linear_cases[experiment]

experiment_data = (slab_length, 0, 0, 1, 0)
simulation_data = (experiment, 60, 2, 1000, 0.05, 'mHNGD')

temperature_data = np.array([[                np.nan,                       0,              slab_length],
                             [                     0, data['low_temperature'], data['high_temperature']],
                             [data['annealing_time'], data['low_temperature'], data['high_temperature']]])
    
initial_hydrogen_data = np.array([[                       0,             slab_length],
                                  [data['initial_hydrogen'], data['initial_hydrogen']]]) 
   
model_parameters = hngd.default.model_parameters

model_parameters['TSSd0'] = 510800
model_parameters['Q_TSSd'] = 45610
model_parameters['TSSp0'] = 66440
model_parameters['Q_TSSp'] = 29630

model_parameters['D0'] = 7.73e-7
model_parameters['Ed'] = 0.47

# Run animation.

anim = hngd.run.animation(temperature_data, initial_hydrogen_data, experiment_data, 
                          simulation_data, model_parameters, interval=2, ylim=300)