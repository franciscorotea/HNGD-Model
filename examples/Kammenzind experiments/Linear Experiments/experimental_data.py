import numpy as np

"""

Experimental data for Kammenzind experiments, as published in Merlino 
J. T., 'Experiments in hydrogen distribution in thermal gradients calculated
using BISON', M. Eng. Thesis, The Pennsylvania State University (2019).

"""

# Linear cases: Appendix A, Experiment 1 Data Sets (pp. 50 to 63).

# 'A' is alpha-annealed Zircaloy-4
# 'B' is Zircaloy-1%Niobium 
# 'C' is Zircaloy-2.5%Niobium

# Specimens were created by cutting larger specimens in half (if both halves 
# were analyzed, 'a' and 'b' corresponds to each one).

linear_cases = {'A26a': {'annealing_time': 27,
                         'low_temperature': 260,
                         'high_temperature': 427,
                         'heat_of_transport': 30214}, 
                'A27': {'annealing_time': 27,
                        'low_temperature': 260,
                        'high_temperature': 427,
                        'heat_of_transport': 28008},
                'A45': {'annealing_time': 77,
                        'low_temperature': 260,
                        'high_temperature': 427,
                        'heat_of_transport': 29185},
                'A46': {'annealing_time': 77,
                        'low_temperature': 260,
                        'high_temperature': 427,
                        'heat_of_transport': 23200},
                'A26b': {'annealing_time': 54,
                         'low_temperature': 260,
                         'high_temperature': 482,
                         'heat_of_transport': 28655},
                'C13': {'annealing_time': 54,
                        'low_temperature': 260,
                        'high_temperature': 482,
                        'heat_of_transport': 11886},
                'A09a': {'annealing_time': 15,
                         'low_temperature': 316,
                         'high_temperature': 482,
                         'heat_of_transport': 28753},
                'A10a': {'annealing_time': 15,
                         'low_temperature': 316,
                         'high_temperature': 482,
                         'heat_of_transport': 30055},
                'A09b': {'annealing_time': 32,
                         'low_temperature': 316,
                         'high_temperature': 482,
                         'heat_of_transport': 33581},
                'A10b': {'annealing_time': 32,
                         'low_temperature': 316,
                         'high_temperature': 482,
                         'heat_of_transport': 30873},
                'A31': {'annealing_time': 26,
                        'low_temperature': 316,
                        'high_temperature': 538,
                        'heat_of_transport': 22878},
                'C14a': {'annealing_time': 26,
                         'low_temperature': 316,
                         'high_temperature': 538,
                         'heat_of_transport': 20824},
                'B14': {'annealing_time': 32,
                        'low_temperature': 316,
                        'high_temperature': 482,
                        'heat_of_transport': 32596},
                'C14b': {'annealing_time': 32,
                         'low_temperature': 316,
                         'high_temperature': 482,
                         'heat_of_transport': 31662},
                'A11a': {'annealing_time': 9,
                         'low_temperature': 371,
                         'high_temperature': 538,
                         'heat_of_transport': 30195},
                'A12a': {'annealing_time': 9,
                         'low_temperature': 371,
                         'high_temperature': 538,
                         'heat_of_transport': 30647},
                'B15': {'annealing_time': 14,
                        'low_temperature': 371,
                        'high_temperature': 538,
                        'heat_of_transport': 28065},
                'C15a': {'annealing_time': 14,
                         'low_temperature': 371,
                         'high_temperature': 538,
                         'heat_of_transport': 25659},
                'A12b': {'annealing_time': 18,
                         'low_temperature': 371,
                         'high_temperature': 538,
                         'heat_of_transport': 27450},
                'C15b': {'annealing_time': 18,
                         'low_temperature': 371,
                         'high_temperature': 538,
                         'heat_of_transport': 28057},
                'A11b': {'annealing_time': 20,
                         'low_temperature': 371,
                         'high_temperature': 593,
                         'heat_of_transport': 15595},
                'C18': {'annealing_time': 20,
                        'low_temperature': 371,
                        'high_temperature': 593,
                        'heat_of_transport': 18547},
                'A13a': {'annealing_time': 6,
                         'low_temperature': 427,
                         'high_temperature': 621,
                         'heat_of_transport': 21773},
                'A14': {'annealing_time': 6,
                        'low_temperature': 427,
                        'high_temperature': 621,
                        'heat_of_transport': 20742},
                'C16a': {'annealing_time': 14,
                         'low_temperature': 427,
                         'high_temperature': 593,
                         'heat_of_transport': 20056},
                'B16': {'annealing_time': 14,
                        'low_temperature': 427,
                        'high_temperature': 593,
                        'heat_of_transport': 27868},
                'A13': {'annealing_time': 12,
                        'low_temperature': 427,
                        'high_temperature': 649,
                        'heat_of_transport': 21773},
                'C16b': {'annealing_time': 12,
                         'low_temperature': 427,
                         'high_temperature': 649,
                         'heat_of_transport': 14716}}


for experiment in linear_cases:

    sample_length, midpoint_location, midpoint_temperature, hydrogen_concentration = np.loadtxt(f'Data/{experiment}.txt', 
                                                                                                usecols=range(1,5), unpack=True)
    
    linear_cases[experiment]['annealing_time'] = linear_cases[experiment]['annealing_time']*24*60 # days to minutes
    linear_cases[experiment]['sample_length'] = sample_length*1e-2 # cm to m
    linear_cases[experiment]['midpoint_location'] = midpoint_location*1e-2 # cm to m
    linear_cases[experiment]['midpoint_temperature'] = midpoint_temperature # °C
    linear_cases[experiment]['hydrogen_concentration'] = hydrogen_concentration # wt.ppm
    linear_cases[experiment]['initial_hydrogen'] = np.dot(hydrogen_concentration, sample_length)/np.sum(sample_length) # wt.ppm