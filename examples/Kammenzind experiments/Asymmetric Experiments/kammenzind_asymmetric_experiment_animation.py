import hngd.run
import hngd.default
import experimental_data
import numpy as np

# Choose an experiment from the data base. The valid options for asymmetric 
# cases are: A54, A56, A53, A55a, A55b.

experiment = 'A56'

# Input parameters.

slab_length = 3.810*1e-2
hydride_rim = 0.100*1e-2

data = experimental_data.asymmetric_cases[experiment]

experiment_data = (slab_length, 0, 0, 1, 0)
simulation_data = (experiment, 80, 2, 50000, 0.05, 'mHNGD')

temperature_data = data['temperature_data']
    

initial_hydrogen_data = np.array([[ 0, slab_length-hydride_rim, slab_length-hydride_rim+1e-12,               slab_length],
                                  [10,                      10,      data['initial_hydrogen'], data['initial_hydrogen']]])
   
model_parameters = hngd.default.model_parameters

model_parameters['TSSd0'] = 510800
model_parameters['Q_TSSd'] = 45610
model_parameters['TSSp0'] = 66440
model_parameters['Q_TSSp'] = 29630

model_parameters['D0'] = 7.73e-7
model_parameters['Ed'] = 0.47

# Run animation.

anim = hngd.run.animation(temperature_data, initial_hydrogen_data, experiment_data, 
                          simulation_data, model_parameters, interval=2)