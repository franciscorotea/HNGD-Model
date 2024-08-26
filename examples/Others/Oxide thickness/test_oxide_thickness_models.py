import numpy as np
import matplotlib.pyplot as plt

# TODO calc_transition es time para COCHISE, temperature para los otros. 
# TODO revisar COCHISE, no funciona!

"""The hydrogen flux entering the cladding is proportional to the oxidation 
rate. This flux can be evaluated from the oxide kinetics equations. The 
oxidation kinetics have been formulated using semi-empirical models, as 
detailed in [1]. These models separate two different kinetic behaviors. At 
first, oxidation weight gain kinetics is governed by a cubic rate law. Then, 
the oxidation rate gradually decreases, until at a given oxide thickness
(1.8-2.0 um for Zircaloy-4), a kinetic transition is observed. At this point, 
often referred to as the oxide transition, the oxidation kinetics return to the 
initial value seen at the start of the corrosion of the bare metal. Henceforth, 
the oxidation kinetics can be approximated with a linear rate law.

Ref.:
    [1] SCDAP/RELAP5/MOD3.1 Code Manual Volume IV: MATPRO -- A Library of 
        Materials Properties for Light-Water-Reactor Accident Analysis, Idaho 
        National Engineering Laboratory. NUREG/CR-6150, EGG-2720, Volume IV.
    [2] IAEA, Waterside Corrosion of Zirconium Alloys in Nuclear Power Plants,
        International Atomic Energy Agency, Vienna TECDOC996, 1998.
    [3] Courty O., Motta A. T., Hales J. D, Modeling and Simulation of Hydrogen 
        Behavior in Zircaloy-4 Fuel Cladding. Journal of Nuclear Materials 452, 
        pp. 311-320 (2014). https://doi.org/10.1016/j.jnucmat.2014.05.013.

"""

dt = 0.1        # [days]

# MATPRO model > typical reactor operation at temperatures of 523 to 673 K.

# A_Kc: pre-transition frequency factor (cubic regime)
# A_Kl: post-transition frequency factor (linear regime)
# Q_R_pre: pre-transition activation energy (Q/R)
# Q_R_post: post-transition activation energy (Q/R)

model_params = {'MATPRO': {'A_Kc': 4.976e9,
                           'Q_R_pre': 15660,
                           'A_Kl': 8.288e7,
                           'Q_R_post': 14080},
                'EPRI': {'A_Kc': 1.78e10,
                         'Q_R_pre': 16250,
                         'A_Kl': 8.04e7,
                         'Q_R_post': 13766},
                'COCHISE': {'A_Kc': 11.4e10,
                            'Q_R_pre': 17171,
                            'A_Kl': 4e11,
                            'Q_R_post': 18391}}

def calc_irradiation_enhancement_factor(model, oxide_coolant_interface_temperature, fast_flux=5e16):
    """
    A function to calculate a parameter that describes enhancement of the 
    cladding oxidation rate in a reactor environment. Typical values are around 
    1.5 and 9 for PWR and BWR, respectively.

    Parameters
    ----------
    model : str
        The model to use. Can be 'MATPRO', 'EPRI' or 'COCHISE'.
    oxide_coolant_interface_temperature : float or ndarray
        Oxide-coolant interface temperature [K].
    fast_flux : float
        Fast neutron flux, greater than 0.8-1.0 MeV [neutrons m-2 s-1].

    Returns
    -------
    irradiation_enhacement_factor : float or ndarray
        Enhancement of oxidation rate due to irradiation.

    """
    
    match model:
        
        case 'MATPRO': # [1]
            #E = 2.95 # ranges between 2.07 and 4.87 (average: 2.95) [2]. Bergas usa 1?
            E = 1
            f = 120.3 * np.exp(-0.007118*oxide_coolant_interface_temperature)

        case 'EPRI': # [2]
            E = 1
            M = 7.46e-19 # m2 s neutr-1
            f = 1 + 3.22*(M*fast_flux)**0.24 # fast flux in neutr m-2 s-1

        case 'COCHISE': # irradiation factor is already included in the model
            E = 1
            f = 1
            
    return E*f

def calc_transition_thickness(model, oxide_metal_interface_temperature):
    """
    A function to calculate the kinetic transition between a cubic rate law and
    a linear rate law, according to different models.

    Parameters
    ----------
    model : str
        The model to use. Can be 'MATPRO', 'EPRI' or 'COCHISE'.
    temperature : float or ndarray
        Temperature [K]

    Returns
    -------
    transition_thickness : float or ndarray
        Transition oxide thickness [um].
    """
    
    match model:
        
        case 'MATPRO': # [1]
            A = 7.749
            B = 790/oxide_metal_interface_temperature
            
        case 'EPRI':
            A = 2.14e7
            B = 5417/oxide_metal_interface_temperature + 0.0117*oxide_metal_interface_temperature
            
        case 'COCHISE':
            A = 8.857e10
            B = -921/oxide_metal_interface_temperature + 0.035*oxide_metal_interface_temperature
            
    transition_thickness = A * np.exp(-B)
    
    return transition_thickness

def calc_oxide_thickness(end_time, model, oxide_metal_interface_temperature, oxide_coolant_interface_temperature):

    t = 0
    d = 0
    
    initial_thickness = 0           # [um]
    
    time = []
    thickness = []
    
    transition_thickness = calc_transition_thickness(model, oxide_metal_interface_temperature)
    
    print(f'Transition thickness using the {model} model for {oxide_metal_interface_temperature-273}°C is {np.round(transition_thickness, 2)} um.')
    
    transition = False
    
    irradiation_enhancement_factor = calc_irradiation_enhancement_factor(model, oxide_coolant_interface_temperature)

    while t < end_time:
        
        # Pre-transition oxide film
        
        if d < transition_thickness:
            
            cubic_constant = model_params[model]['A_Kc'] * np.exp(-model_params[model]['Q_R_pre']/oxide_metal_interface_temperature)
            
            # In the EPRI model, irradiation factor is applied to the post-transition region but not to the pre-transition region.
                
            if model == 'MATPRO':
                d = (cubic_constant*irradiation_enhancement_factor*t + initial_thickness**3)**(1/3)
            elif model == 'EPRI':
                d = (cubic_constant*t + initial_thickness**3)**(1/3)
                
        # Post-transition oxide film
        
        else:
            
            if not transition:
                transition_time = t
                transition = True
                print(f'Transition time using the {model} model is = {np.round(transition_time, 2)} days.')
            
            if initial_thickness < transition_thickness:
                linear_constant = model_params[model]['A_Kl'] * np.exp(-model_params[model]['Q_R_post']/oxide_metal_interface_temperature)
                d = transition_thickness + irradiation_enhancement_factor * linear_constant * (t-transition_time)
            else:
                linear_constant = model_params[model]['A_Kl'] * np.exp(-model_params[model]['Q_R_post']/oxide_metal_interface_temperature)
                d = initial_thickness + irradiation_enhancement_factor * linear_constant * t
            
        time.append(t)
        thickness.append(d)
        
        t = t + dt
        
    return time, thickness

fig, ax = plt.subplots()

for model in ('MATPRO', 'EPRI', 'COCHISE'):
    
    time, thickness = calc_oxide_thickness(1500, model, 300+273, 300+273)

    ax.plot(time, thickness)
    
    ax.set_xlabel('Time [d]')
    ax.set_ylabel('Thickness [um]')

plt.show()