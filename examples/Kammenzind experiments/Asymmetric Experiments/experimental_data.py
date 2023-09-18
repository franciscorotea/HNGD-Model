import numpy as np

"""

Experimental data for Kammenzind experiments, as published in Merlino 
J. T., 'Experiments in hydrogen distribution in thermal gradients calculated
using BISON', M. Eng. Thesis, The Pennsylvania State University (2019).

"""

# Asymmetric cases: Appendix A, Experiment 2 Data Sets (pp. 63 to 70).

# All experiments are alpha-annealed Zircaloy-4.

# Specimens were created by cutting larger specimens in half (if both halves 
# were analyzed, 'a' and 'b' corresponds to each one).

asymmetric_cases = {'A54': {'annealing_time': 194,
                            'temperature_data': np.array([[  np.nan, 0*1e-2, 0.254*1e-2, 0.508*1e-2, 0.762*1e-2, 1.334*1e-2, 1.905*1e-2, 2.540*1e-2, 3.175*1e-2, 3.493*1e-2, 3.810*1e-2],
                                                         [        0,  261.7,      265.9,      269.7,      280.7,      319.1,      358.4,      374.0,      349.6,      324.3,      303.1],
                                                         [194*24*60,  261.7,      265.9,      269.7,      280.7,      319.1,      358.4,      374.0,      349.6,      324.3,      303.1]]),
                            'hydride_rim_concentration': 2675,
                            'heat_of_transport': 24767}, 
                    'A56': {'annealing_time': 95,
                            'temperature_data': np.array([[  np.nan, 0*1e-2, 0.254*1e-2, 0.508*1e-2, 0.762*1e-2, 1.334*1e-2, 1.905*1e-2, 2.540*1e-2, 3.175*1e-2, 3.493*1e-2, 3.810*1e-2],
                                                         [       0,   261.1,      264.7,      268.1,      279.1,      314.6,      352.3,      374.9,      350.1,      326.9,      303.3],
                                                         [95*24*60,   261.1,      264.7,      268.1,      279.1,      314.6,      352.3,      374.9,      350.1,      326.9,      303.3]]),
                            'hydride_rim_concentration': 3478,
                            'heat_of_transport': 26566},
                    'A53': {'annealing_time': 150,
                            'temperature_data': np.array([[  np.nan, 0*1e-2, 0.254*1e-2, 0.508*1e-2, 0.762*1e-2, 1.334*1e-2, 1.905*1e-2, 2.540*1e-2, 3.175*1e-2, 3.493*1e-2, 3.810*1e-2],
                                                         [        0,  260.9,      267.9,      276.6,      286.7,      320.1,      360.9,      371.4,      337.5,      305.7,      289.1],
                                                         [150*24*60,  260.9,      267.9,      276.6,      286.7,      320.1,      360.9,      371.4,      337.5,      305.7,      289.1]]),
                            'hydride_rim_concentration': 3542,
                            'heat_of_transport': 32601},
                    'A55a': {'annealing_time': 209,
                             'temperature_data': np.array([[  np.nan, 0*1e-2, 0.254*1e-2, 0.508*1e-2, 0.762*1e-2, 1.334*1e-2, 1.905*1e-2, 2.540*1e-2, 3.175*1e-2, 3.493*1e-2, 3.810*1e-2],
                                                          [        0,  317.2,      319.7,      321.5,      326.7,      342.7,      364.8,      372.9,      364.5,      354.8,      345.5],
                                                          [209*24*60,  317.2,      319.7,      321.5,      326.7,      342.7,      364.8,      372.9,      364.5,      354.8,      345.5]]),
                             'hydride_rim_concentration': 3077,
                             'heat_of_transport': 42107},
                    'A55b': {'annealing_time': 209,
                             'temperature_data': np.array([[  np.nan, 0*1e-2, 0.254*1e-2, 0.508*1e-2, 0.762*1e-2, 1.334*1e-2, 1.905*1e-2, 2.540*1e-2, 3.175*1e-2, 3.493*1e-2, 3.810*1e-2],
                                                          [        0,  317.2,      319.7,      321.5,      326.7,      342.7,      364.8,      372.9,      364.5,      354.8,      345.5],
                                                          [209*24*60,  317.2,      319.7,      321.5,      326.7,      342.7,      364.8,      372.9,      364.5,      354.8,      345.5]]),
                             'hydride_rim_concentration': 3125,
                             'heat_of_transport': 27163}}

for experiment in asymmetric_cases:

    sample_length, midpoint_location, midpoint_temperature, hydrogen_concentration = np.loadtxt(f'Data/{experiment}.txt', 
                                                                                                usecols=range(1,5), unpack=True)
    
    asymmetric_cases[experiment]['annealing_time'] = asymmetric_cases[experiment]['annealing_time']*24*60 # days to minutes
    asymmetric_cases[experiment]['sample_length'] = sample_length*1e-2 # cm to m
    asymmetric_cases[experiment]['midpoint_location'] = midpoint_location*1e-2 # cm to m
    asymmetric_cases[experiment]['midpoint_temperature'] = midpoint_temperature # °C
    asymmetric_cases[experiment]['hydrogen_concentration'] = hydrogen_concentration # wt.ppm
    asymmetric_cases[experiment]['initial_hydrogen'] = asymmetric_cases[experiment]['hydride_rim_concentration'] # wt.ppm