import hngd.run
import numpy as np

"""
A sample script to animate a thermal history.
"""

temperature_data = np.array([[np.nan,   0, 0.025/2, 0.025],
                             [     0,  20,      50,    20],
                             [    38, 400,     200,   400],
                             [    48, 400,     400,   400],
                             [    68, 400,     400,   400],
                             [   108,  20,      80,    20]])

# Run simulation

anim = hngd.run.thermal_history_animation(temperature_data, n_nodes=10, dt=0.6,
                                          interval=50)