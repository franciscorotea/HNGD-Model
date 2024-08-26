"""A module with an implementation of the Hydride Nucleation-Growth-Dissolution
(HNGD) model, that simulates hydrogen behaviour in Zirconium alloys in 1D 
geometry using an explicit finite difference scheme (Euler method).

Ref.:
    
    [1] Lacroix E., 'Modeling zirconium hydride precipitation and dissolution
        in zirconium alloys', PhD Thesis, The Pennsylvania State University 
        (2019).
    [2] Passelaigue F., 'Hydride Nucleation-Growth-Dissolution model: 
        Implementation in BISON', MSc. Thesis, The Pennsylvania State 
        University (2020).
    [3] Passelaigue F., Lacroix E., Pastore G. and Motta A. T., 'Implementation 
        and validation of the Hydride Nucleation-Growth-Dissolution (HNGD) 
        model in BISON', Journal of Nuclear Materials 544, pp. 152683 (2021).
        https://doi.org/10.1016/j.jnucmat.2020.152683
    [4] Passelaigue F., Simon, P. A. and Motta, A. T., 'Predicting the hydride
        rim by improving the solubility limits in the Hydride Nucleation-Growth
        -Dissolution (HNGD) model', Journal of Nuclear Materials 558, 
        pp. 153363 (2022).
        https://doi.org/10.1016/j.jnucmat.2021.153363
    [5] Bin Seo S., Duchnowski E. M., Motta A. T., Kammenzind B. F. and Brown 
        N. R.,'Sensitivity analysis for characterizing the impact of HNGD model
        on the prediction of hydrogen redistribution in Zircaloy cladding using
        BISON code', Nuclear Engineering and Design 393, pp. 111813 (2022).
        https://doi.org/10.1016/j.nucengdes.2022.111813
    [6] Lee C. and Lee, Y., 'Simulation of hydrogen diffusion along the axial 
        direction in zirconium cladding tube during dry storage', Journal of
        Nuclear Materials 579, pp. 154352 (2023).
        https://doi.org/10.1016/j.jnucmat.2023.154352
        
"""

import numpy as np
from numba import njit
#from numba import config
#config.DISABLE_JIT=True

###############################################################################

# TODO > 

# 1- Incluir hydrostatic_stress (solo habria que hacer un input analogo a 
# temperature_data, su efecto ya esta incluido en la funcion 
# calc_diffusion_flux). Revisar unidades (Pa?).

# 2- Se podria cambiar la discretizacion espacial para que se pueda hacer mas 
# fina en los extremos. Por ejemplo:
    
# =============================================================================
# def denseboundspace(size=30, start=0, end=10, a=2, b=5):
#     x = np.linspace(0, 1, size)
#     return start + beta.cdf(x, a, b) * (end-start)
# 
# n_x = denseboundspace()
# 
# plt.plot(n_x, np.zeros_like(n_x), 'o')
# 
# plt.show()
# =============================================================================

# donde 'size' es el numero de elementos, 'start' y 'end' son los valores entre
# los cuales se quiere hacer la discretizacion (0 y slab_length), y 'a', 'b' 
# controlan que tan "densa" es la discretizacion en los extremos izquierdo y 
# derecho, respectivamente. La funcion 'beta.cdf' esta en la libreria 
# 'scipy.stats' (es la distribucion beta). Entonces 'dx' pasaria a ser un array.

# 3- Poner el modelo mHNGD en una funcion (refactor > hay cosas que se repiten, 
# como C_delta y lever_rule_Cp que se calculan en precipitation_dissolution).

# 4- En las animaciones no funciona repeat=True. Seguramente tiene que ver con 
# el pass by reference que usa Python. Por ejemplo, en 
# 'Examples/animation_example.py', si se usa repeat=True solamente funciona 
# bien la repeticion de TSSp, el resto solo funciona bien en la primera pasada.
# Pero cambiando en las ultimas lineas del if model = 'mHNGD' de 'simulate' 
# (lineas 1010 y 1011) lo siguiente:
    
#    TSSd = np.copy(TSSd_eff)
#    TSSp = np.copy(TSSp_eff)

# Pasa a funcionar bien el TSSd. Habria que hacer algo asi para Css, Cp, 
# temperature_profile, etc etc.

# 6- En 'calc_kinetic_parameters', revisar por que es necesario limitar a uno 
# el parametro de nucleacion. Sale del codigo original de Lacroix, para evitar 
# la sobreprecipitacion? Probe sacando esa condicion y no modifica los 
# resultados, pero hace muy lentos los calculos (K_D ~ 2000, 1/K_D = dt = 0.0005).

###############################################################################

# Physics constants

k_B = 8.617333262e-5 # Boltzmann constant [eV K-1]
R = 8.31446261815324 # Gas constant [J mol-1 K-1]
m_H = 1.008 # Molar mass of Hydrogen [g/mol]
m_Zr = 91.224 # Molar mass of Zirconium [g/mol]

@njit
def profile_changed(current_profile, reference_profile, criterion):
    """
    Returns True if the current profile is significantly different to the
    saved reference profile, according to a criterion fixed by the user. Used 
    to see if the hydrogen or temperature profile changed significantly, in 
    order to discern if the output will be saved or not. If True, the current 
    profile would take over as reference.
    
    For instance, criterion = 0.05 represents a 5% change between the current
    profile and the reference profile.

    Parameters
    ----------
    current_profile : 1d numpy array 
        Current hydrogen/temperature profile.
    reference_profile : 1d numpy array 
        Saved hydrogen/temperature profile.
    criterion : float greater than 0
        Criterion to compare current and reference profiles. Typical values
        are 0.01 or 0.05 (1% or 5% change).

    Returns
    -------
    profile_changed : bool
        Returns True if the profile changed according to the set criterion, and
        False otherwise.

    """
    
    relative_difference = np.abs((current_profile-reference_profile)/reference_profile)
    average_relative_difference = np.mean(relative_difference)
    
    return average_relative_difference > criterion

@njit
def correct_input_positions(input_positions, x):
    """
    Correct temperature/hydrogen input positions so that each position given
    by the user is equal to a node position. This is done by selecting the 
    closest values to the input position array in the spatial mesh. Given a 
    sufficiently fine mesh, the approximated value should not be too far to the 
    one given by the user.

    Parameters
    ----------
    input_positions : 1d numpy array
        The temperature/hydrogen input position array given by the user [m].
    x : 1d numpy array
        Spatial mesh [m].

    Returns
    -------
    corrected_input_positions : 1d numpy array
        The corrected temperature/hydrogen input position array, with positions
        that belong to the mesh x [m].

    """
    
    return np.array([find_nearest(position, x) for position in input_positions])

@njit
def find_nearest(value, array):
    """
    Finds the nearest value in an array. For instance, if the array is:
        
        x = np.array([1.2, 3, 5.6, 10])
        
    And value = 5, the function would return 5.6.

    Parameters
    ----------
    value : float
        The value you want to find the closest value in the array.
    array : 1d numpy array
        A one dimensional numpy array.

    Returns
    -------
    closest_value : float
        Closest value to `value` in the array.

    """
    
    idx = (np.abs(array - value)).argmin()
    
    return array[idx]

@njit
def time_interpolation(t, input_time_stamps, input_temperatures):
    """
    Performs linear interpolation to find the temperatures at the positions
    indicated in the input at time t.

    Parameters
    ----------
    input_time_stamps : 1d numpy array of size `n`
        Array with the time stamps of the input temperature history [s].
    input_temperatures : 2d numpy array of size `(n, k)`
        Array with the input temperature history of the sample [K]. For each
        of the `n` time stamps, it has `k` values representing the known 
        temperatures at different positions of the sample.
    t : float
        The time at which you want to find the corresponding temperature 
        stamps [s].

    Returns
    -------
    interpolated_temperatures : 1d numpy array of size `k`
        The temperatures at each input position of the sample, interpolated 
        for the time `t`.

    """
    
    n_positions = np.shape(input_temperatures)[1]
    
    return np.array([np.interp(t, input_time_stamps, input_temperatures[:,i]) for i in range(n_positions)])

@njit
def spatial_interpolation(x, input_positions, values):
    """
    Calculates the temperature/hydrogen profile at each node of the mesh using 
    linear interpolation from the values known at the input positions.

    Parameters
    ----------
    x : 1d numpy array of size `q`
        Array with the positions determined by the mesh [m].
    input_positions : 1d numpy array of size `k`
        Array with the positions determined by the input temperature/hydrogen 
        history.
    values : 1d numpy array of size `k`
        Array with the values of temperature/hydrogen at each input position.

    Returns
    -------
    interpolated_profile : 1d numpy array of size `q`
        Array with the interpolated values of temperature/hydrogen at each node
        of the mesh.

    """
    
    return np.interp(x, input_positions, values)

@njit
def calc_diffusion_coefficient(D0, Ed, temperature):
    """
    Calculates the diffusion coefficient at each node of the sample.

    Parameters
    ----------
    D0 : float
        Pre-exponential term of the diffusion coefficient [m2 s-1].
    Ed : float
        Activation energy of the diffusion coefficient [eV at-1]
    temperature : 1d numpy array
        Temperature profile [K].

    Returns
    -------
    diffusion_coefficient : 1d numpy array
        Diffusion coefficient at each node of the sample [m2 s-1].

    """
    
    return D0 * np.exp(-Ed/(k_B*temperature))

@njit
def calc_terminal_solid_solubility(TSSp0, Q_TSSp, TSSd0, Q_TSSd, temperature,
                                   x=None, liner_width=0.0, 
                                   liner_solubility_factor=1.0):
    """
    Calculates terminal solid solubility for precipitation (TSSp) and 
    dissolution (TSSd), also called solubility and supersolubility limits. 
    
    A liner (a thin, protective layer that is applied to the inner surface of 
    the cladding) can be implemented through a local difference in solubility: 
    liners are usually made of pure Zr and don't have much oxygen, causing the 
    solubility in the liner to be lower than in the rest of the cladding. For
    a more comprehensive description, see Section 5.3.2: Segregation of 
    Hydrides in Liners in Ref. [1].
    
    Parameters
    ----------
    TSSp0 : float
        Pre-exponential term for TSSp [wt ppm].
    Q_TSSp : float
        Activation energy for TSSp [J mol-1].
    TSSd0 : float
        Pre-exponential term for TSSd [wt pmm]
    Q_TSSd : float
        Activation energy for TSSd [J mol-1].
    temperature : 1d numpy array
        Temperature profile [K].
    x :  1d numpy array, optional
        Array with the positions determined by the mesh [m].
    liner_width : float, optional
        The width of the liner [m].
    liner_solubility_factor : float between 0 and 1, optional
        Scalar to multiply the solubility to account for a liner. It was
        observed that 0.95-0.98 gives good results.

    Returns
    -------
    TSSp : 1d numpy array
        Terminal solid solubility for precipitation at each node [wt ppm].
    TSSd : 1d numpy array
        Terminal solid solubility for dissolution at each node [wt ppm].

    """
    
    TSSp = TSSp0 * np.exp(-Q_TSSp/(R*temperature))
    TSSd = TSSd0 * np.exp(-Q_TSSd/(R*temperature))
    
    if x is not None and liner_width > 0:
        
        liner = x < liner_width
        
        TSSp[liner] = TSSp0 * np.exp(-Q_TSSp/(R*temperature[liner])) * liner_solubility_factor
        TSSd[liner] = TSSd0 * np.exp(-Q_TSSd/(R*temperature[liner])) * liner_solubility_factor   
    
    return TSSp, TSSd

@njit
def set_initial_hydrogen_distribution(n_nodes, TSSd, initial_hydrogen_profile):
    """
    Set the initial hydrogen distribution across each node of the sample.
    Steady state is assumed to determine how much hydrogen content is present
    in solid solution and how much is present as hydrides.

    Parameters
    ----------
    n_nodes : int
        Number of nodes in the sample [/].
    TSSd : 1d numpy array of size `q`
        Terminal solid solubility for dissolution at each node [wt ppm].
    initial_hydrogen_profile : 1d numpy array of size `q`
        The initial profile of total hydrogen at each node of the sample 
        [wt ppm].

    Returns
    -------
    Css : 2d numpy array of size `(q,2)`
        Concentration of hydrogen in solid solution [wt ppm].
    Cp : 1d numpy array of size `q`
        Concentration of hydrogen in hydrides [wt ppm].

    """
    
    # Arrays for hydrogen in solid solution, Css [wt.ppm] and hydrogen in 
    # hydrides, Cp [wt.ppm].
    
    Css = np.zeros((n_nodes+3, 2))
    Cp = np.zeros(n_nodes+3)
    
    # Create a mask to compare the total hydrogen initial profile with TSSd.
    
    mask = TSSd < initial_hydrogen_profile
    
    # At the start of the experiment, steady state is assumed. Then hydrogen in
    # solid solution equals the solvus TSSd, and the rest is present as 
    # hydrides.
    
    Css[mask, 0] = TSSd[mask]
    Cp[mask] = initial_hydrogen_profile[mask] - TSSd[mask]
    
    # If there is less hydrogen than the solvus TSSd, all hydrogen is in solid 
    # solution, and there are no hydrides present.
    
    Css[~mask, 0] = initial_hydrogen_profile[~mask]
    Cp[~mask] = 0 
    
    return Css, Cp

@njit
def ppm_to_atomic_fraction(x):
    """
    Converts hydrogen concentration `x` from parts per million (ppm) to atomic 
    fraction.

    Parameters
    ----------
    x : 1d numpy array
        Hydrogen concentration [wt ppm].

    Returns
    -------
    y : 1d numpy array
        Hydrogen concentration [atomic fraction].

    """

    return x / (m_H * (x/m_H + (1e6-x)/m_Zr))

@njit
def calc_kinetic_parameters(Eth0, Eth1, Eth2, Eth3, K_G_mob0, K_G_th0, K_D_0, 
                            K_N_0, Ed, Eg, temperature, Cp, Css, TSSd):
    """
    Calculation of kinetic parameters for nucleation (N), growth (G) and 
    dissolution (D).
    
    Parameters
    ----------
    Eth0, Eth1, Eth2, Eth3 : floats
        Coefficients of the 3rd degree polynomial used to express the formation 
        energy of delta-hydrides
    K_G_mob0, K_G_th0, K_D_0, K_N_0 : floats
        Pre-exponential factors for growth (diffusion and reaction controlled), 
        nucleation and dissolution kinetics [s-1].
    Ed, Eg : floats
        Activation energy of the diffusion coefficient and diffusion driven 
        growth [eV at-1].
    temperature : 1d numpy array of size `q`
        Temperature profile [K].
    Cp : 1d numpy array of size `q`
        Concentration of hydrogen in hydrides [wt ppm].
    Css : 1d numpy array of size `q`
        Concentration of hydrogen in solid solution [wt ppm].
    TSSd : 1d numpy array of size `q`
        Terminal solid solubility for dissolution at each node [wt ppm].

    Returns
    -------
    K_N, K_G, K_D : 1d numpy array of size `q`
        Kinetic parameters for nucleation, growth, and dissolution [s-1].
    max_K_N, max_K_G, max_K_D: floats
        Maximum rate of nucleation, growth and dissolution [s-1].

    """
    
    # Calculate the formation energy of delta hydrides at a given temperature
    # [eV/atom]. Used for the nucleation and growth of hydrides.
    
    Eth = -Eth0 + Eth1*temperature - Eth2*temperature**2 + Eth3*temperature**3
    
    ############################### LEVER RULE ################################
    
    # Convert mass fraction of total hydrogen (Co=Cp+Css), hydrides (Cp) and
    # TSSd from ppm to atomic fraction of sample studied.
    
    TSSd_at = ppm_to_atomic_fraction(TSSd)
    Co_at = ppm_to_atomic_fraction(Cp + Css)
    Cp_at = ppm_to_atomic_fraction(Cp)
    
    # Calculating the alpha+delta/delta boundary value in the phase diagram.
    
    C_delta = -9.93e-11*temperature**3 + 8.48e-8*temperature**2 - 5.73e-5*temperature + 0.623
    
    # Calculating the lever rule of the phase diagram
    
    lever_rule_Co = (Co_at - TSSd_at) / (C_delta - TSSd_at)
    lever_rule_Cp = Cp_at / (C_delta - TSSd_at)
    
    ################################# GROWTH ##################################
    
    # Defining the pre-exponential factor of growth kinetics to be linearly 
    # dependent on hydrogen content and hydride content.
        
    Ko = K_G_mob0 * lever_rule_Co * (1-lever_rule_Cp) # Diffusion-controlled pre-exp
    K1 = K_G_th0 * lever_rule_Co * (1-lever_rule_Cp) # Reaction-controlled pre-exp
        
    K_mob = Ko * np.exp(-Eg/(k_B*temperature)) # Diffusion-controlled kinetic factor
    K_ther = K1 * np.exp(-Eth/(k_B*temperature)) # Reaction-controlled kinetic factor
    
    # The slowest process is the process controlling the reaction > parallel 
    # process averaging.
            
    K_G = 1 / (1/K_mob + 1/K_ther)
    
    ############################### DISSOLUTION ###############################
    
    # Defining dissolution kinetics parameter: diffusion-controlled.

    K_D = K_D_0 * np.exp(-Ed/(k_B*temperature))
    
    ################################ NUCLEATION ###############################

    # Making the maximum nucleation rate 1, to avoid having overprecipitation.
    
    Kn0 = K_N_0 * (1 - lever_rule_Cp)
    
    K_N = np.minimum(np.ones_like(Cp), 
                     (17000-Cp)/17000*Kn0*np.exp(-Eth/(k_B*temperature)))
    
    ############################# MAXIMUM VALUES ##############################
    
    # Calculating the maximum growth rate, dissolution rate and nucleation rate
    # in the sample. Used to define an adequate timestep.
    
    max_K_N = np.max(np.maximum(K_N, np.zeros_like(K_N)))
    max_K_G = np.max(np.maximum(K_G, np.zeros_like(K_G)))
    max_K_D = np.max(np.maximum(K_D, np.zeros_like(K_D)))
    
    if max_K_N == 0:
        max_K_N = 1
        
    if max_K_G == 0:
        max_K_G = 1
    
    if max_K_D == 0:
        max_K_D = 1
    
    return K_N, K_G, K_D, max_K_N, max_K_G, max_K_D

@njit
def calc_diffusion_flux(n_nodes, flux_left, flux_right, Q_star, V_star, dx, temperature, D, 
                        Css, hydrostatic_stress):
    """
    Calculates the hydrogen atoms diffusion flux. Atom migration is driven by
    concentration gradients (Fick's diffusion), temperature gradients (Soret 
    effect) and stress gradients. 
    
    Fick's law redistributes hydrogen from areas of high-concentration to areas
    of low-concentration (tends to homogenize the distribution), the Soret 
    effect redistributes hydrogen from high-temperature regions to 
    low-temperature regions, and stress gradients redistributes hydrogen from 
    areas of low tensile stress to high tensile stress.

    Parameters
    ----------
    n_nodes : int
        Number of nodes in the sample [/].
    flux_left : float
        Boundary condition: flux at the left side of the sample [wt.ppm m-2 s-1].
    flux_right : float
        Boundary condition: flux at the right side of the sample [wt.ppm m-2 s-1].
    Q_star : float
        Heat of transport of hydrogen in zirconium [J mol-1].
    V_star : float
        Partial molar volume of hydrogen in the zirconium matrix [m3 mol-1].
    dx : float
        Size of the nodes in the sample [m].
    temperature : 1d numpy array of size `q`
        Temperature profile [K].
    D : 1d numpy array of size `q`
        Diffusion coefficient at each node of the sample [m2 s-1].
    Css : 1d numpy array of size `q`
        Concentration of hydrogen in solid solution [wt ppm].
    hydrostatic_stress : 1d numpy array of size `q`
        Applied hydrostatic stress [Pa].

    Returns
    -------
    diffusion_flux : 1d numpy array of size `q-1`
        Diffusion flux [wt.ppm m-2 s-1].

    """
    
    diffusion_flux = np.zeros(n_nodes+2) # Hydrogen diffusion flux
    
    # It would be easier to use diff_u = np.diff(u[1:-1,1]), but it wouldn't
    # work with numba because it does an internal reshape which is only 
    # supported for contiguous arrays.
    
    d_Css = Css[2:-1] - Css[1:-2]
    d_temperature = temperature[2:-1] - temperature[1:-2]
    d_stress = hydrostatic_stress[2:-1] - hydrostatic_stress[1:-2]
    
    diffusion_flux[1:-1] = - D[1:-2] * d_Css / dx \
                           - D[1:-2] * Css[1:-2] * Q_star / (R*temperature[1:-2]**2) * d_temperature / dx \
                           + D[1:-2] * Css[1:-2] * V_star / (R*temperature[1:-2]) * d_stress / dx
    
    # Boundary conditions of hydrogen flux
    
    diffusion_flux[0] = flux_left
    diffusion_flux[-1] = flux_right
    
    return diffusion_flux

@njit
def calc_reaction_type(Css, Cp, TSSd, TSSp, activate_nucleation, 
                       activate_dissolution, activate_growth):
    
    # Determine reaction type
    
    nucleation = False
    growth = False
    dissolution = False
    nucleation_and_growth = False
    
    if np.any(Css>TSSd):
        if np.any(Css>TSSp) and activate_nucleation:
            nucleation = True
        if np.any(Cp>0) and activate_growth:
            growth = True
        if np.any(np.logical_and(Css>TSSp, Cp>0)) and activate_nucleation and activate_growth:
            nucleation_and_growth = True
    elif np.any(Css<TSSd):
        if np.any(Cp>0) and activate_dissolution:
            dissolution = True
    
    reaction_type = np.array([True, nucleation, growth, dissolution, nucleation_and_growth, True])
    
    return reaction_type

@njit
def time_step(dx, max_K_G, max_K_N, max_K_D, Q_star, dt_div, 
              diffusion_coefficient, temperature, test_t_set, reaction_type):
    
    """Determine the stable time step from the forward Euler criterion. The
    applied time step is adaptively changed, as temperature and concentration
    change with time.
    
    Every node might undergo different kinetic processes (i. e. growth in one 
    node and dissolution in another node), requiring different stable time 
    steps. Then, it is chosen the minimum required time step of all positions.

    """
    
    D = diffusion_coefficient[1:-2]
    T = temperature[1:-2]
    dT = temperature[2:-1] - temperature[1:-2]
    
    dt_flux = np.min(0.5 * dx**2 / (D*(1+Q_star*dT/(2*R*T**2))))
    dt_nucleation = 1/max_K_N
    dt_growth = 1/(5*max_K_G)
    dt_dissolution = 1/max_K_D
    dt_nuc_and_growth = 1/(max_K_N+5*max_K_G)
    
    if test_t_set == 0:
        dt_test_t_set = 1e12
    else:
        dt_test_t_set = test_t_set
    
    # Compute the minimum time step that ensures a stable solution.
    
    dt_array = np.array([dt_flux, dt_nucleation, dt_growth, dt_dissolution, dt_nuc_and_growth, dt_test_t_set])
    dt_array = dt_array[reaction_type]
    dt = np.min(dt_array/dt_div)
    
    return dt

@njit
def precipitation_dissolution(n_nodes, dt, Css, Cp, TSSd, TSSp, K_N, K_G, K_D, avrami_parameter,
                              activate_nucleation, activate_growth, activate_dissolution):
    """
    Implementation of the hydride precipitation (nucleation+growth) and 
    dissolution mechanisms.

    Parameters
    ----------
    n_nodes : int
        Number of nodes in the sample [/].
    dt : float
        Time step [s].
    Css : 2d numpy array of size `(q,2)`
        Concentration of hydrogen in solid solution [wt ppm].
    Cp : 1d numpy array of size `q`
        Concentration of hydrogen in hydrides [wt ppm].
    TSSd : 1d numpy array of size `q`
        Terminal solid solubility for dissolution at each node [wt ppm].
    TSSp : 1d numpy array of size `q`
        Terminal solid solubility for precipitation at each node [wt ppm].
    K_N, K_G, K_D : 1d numpy array of size `q`
        Kinetic parameters for nucleation, growth, and dissolution [s-1].
    avrami_parameter : float
        Avrami parameter for platelets [/].
    activate_nucleation, activate_growth, activate_dissolution : bool
        Activate nucleation, growth and/or dissolution kinetics. Only used
        for verification purposes.

    Returns
    -------
    Css : 2d numpy array of size `(q,2)`
        Concentration of hydrogen in solid solution [wt ppm].
    Cp : 1d numpy array of size `q`
        Concentration of hydrogen in hydrides [wt ppm].
    growth : float
        Proportion of precipitation due to growth.
    nucleation : float
        Proportion of precipitation due to nucleation.

    """
    
    # These variables calculate what portion of the precipitation was from 
    # nucleation and what portion from growth.
    
    nucleation = 0
    growth = 0
    
    # These variables keep track of the kinetics of precipitation and
    # dissolution of hydrides.
    
    nucleation_rate = 0
    growth_rate = 0
    dissolution_rate = 0
    
    for i in range(n_nodes+3):
        
        ####################### NUCLEATION and GROWTH #########################
    
        if Css[i,1] > TSSd[i]:
            
            # Nucleation
            
            if Css[i,1] > TSSp[i] and activate_nucleation:
                nucleation_rate = K_N[i] * (Css[i,1]-TSSp[i]) # dCss/dt [ppm/s]
            else:
                nucleation_rate = 0
            
            # Growth
            
            if Cp[i] > 0 and activate_growth:
            
                # Calculate alpha, the advancement of the reaction 
                # (1-x)H + xZr -> Zr_xH_(1-x). It is achieved (and = 1) when
                # Css = TSSd, i.e. all the hydrogen content above TSSd is in
                # the form of hydrides.
                
                alpha = Cp[i] / (Cp[i] + Css[i,1] - TSSd[i])
                    
                if alpha >= 1:
                    eps = np.finfo(np.float64).eps
                    alpha = alpha = 1 - eps
                
                # Calculate growth rate.
                            
                growth_rate = K_G[i] * (Cp[i] + Css[i,1] - TSSd[i]) * avrami_parameter * (1-alpha) * (-np.log(1-alpha))**(1-1/avrami_parameter) # dCss/dt [ppm/s]
                
            # Calculating hydride content and hydrogen content in solid 
            # solution in each node.
            
            Cp[i] = Cp[i] + (nucleation_rate+growth_rate)*dt
            Css[i,1] = Css[i,1] - (nucleation_rate+growth_rate)*dt
            
            nucleation = nucleation + nucleation_rate*dt/(n_nodes+3)
            growth = growth + growth_rate*dt/(n_nodes+3)
            
        ############################## DISSOLUTION ############################    
        
        elif Css[i,1] < TSSd[i]:
            
            # Hydrides are present.
            
            if Cp[i] > 0 and activate_dissolution:
                
                dissolution_rate = K_D[i] * (TSSd[i] - Css[i,1]) # dCss/dt [ppm/s]
                
                if dissolution_rate > 0:
                    
                    Css[i,1] = Css[i,1] + dissolution_rate*dt
                    Cp[i] = Cp[i] - dissolution_rate*dt
                    
                    # Failsafe for negative hydride concentration.
                    
                    if Cp[i] < 0:
                        Css[i,1] = Css[i,1] + Cp[i]
                        Cp[i] = 0
                
        ######## NO PRECIPITATION, NO DISSOLUTION (PHASE EQUILIBRIUM) #########
        
        else:
            
            Css[i,1] = Css[i,1]
            Cp[i] = Cp[i]            
            
    return Css, Cp, growth, nucleation

@njit
def run_thermal_history(temperature_data, n_nodes, dt):
    """A simple generator used to animate the thermal history of the sample."""
    
    input_temperature_positions = temperature_data[0, 1:]
    input_time_stamps = temperature_data[1:, 0]
    input_temperature = temperature_data[1:, 1:]
    
    t_max = input_time_stamps[-1]
    
    # Calculate spatial discretization
    
    slab_length = input_temperature_positions[-1]
    dx = slab_length/n_nodes
    x = np.linspace(-dx, slab_length+dx, n_nodes+3)
    
    # Correct input positions so that each position is equal to a node position
    
    input_temperature_positions = correct_input_positions(input_temperature_positions, x)
    
    t = 0
    
    while t < t_max:
        
        # Use linear interpolation to find the temperature values at the input 
        # positions in the sample for time t.
        
        temperatures_at_input_pos = time_interpolation(t, input_time_stamps, input_temperature)
        
        # Compute the temperature profile in the sample
        
        temperature_profile = spatial_interpolation(x, input_temperature_positions, temperatures_at_input_pos)
        
        yield t, x, temperature_profile, input_temperature_positions, temperatures_at_input_pos
        
        t = t + dt

@njit
def simulate(temperature_data, hydrogen_data, experiment_data, simulation_data, model_parameters,
             activate_nucleation=True, activate_growth=True, activate_dissolution=True):
    """
    Main generator to perform a simulation.

    Parameters
    ----------
    temperature_data : 2d numpy array of size `(n+1, q+1)`
        Thermal history of the sample. The first row has the `q` positions [m] 
        and the first column has the `n` time stamps [minutes]. Temperature 
        in [°C].
    hydrogen_data : 2d numpy array of size `(2, s)`
        Hydrogen profile at the start of the experiment. The first row has the 
        positions at which the hydrogen concentration is defined [m]. The 
        second row has the hydrogen concentration at each of these 
        positions [wt.ppm].
    experiment_data : 5-tuple (float, float, float, float, float)
        Data from the experiment: slab length [m], flux at the left of the 
        sample [?], flux at the right of the sample [?], liner solubility
        factor [/] and liner width [m].
    simulation_data : 5-tuple (str, int, float, float, float, str)
        Data from the simulation: name, number of nodes [/], dt_div [s], 
        test_t_set [s], output_criterion [/] and model ['HNGD' or 'mHNGD'].
    model_parameters : dict
        A dictionary with model parameters for the HNGD model.
    activate_nucleation, activate_growth, activate_dissolution : bool
        Activate nucleation, growth and/or dissolution kinetics. Only used
        for verification purposes (defaults to True).

    Yields
    ------
    Css : 1d numpy array of size `q`
        Concentration of hydrogen in solid solution [wt ppm].
    Cp : 1d numpy array of size `q`
        Concentration of hydrogen in hydrides [wt ppm].
    TSSp : 1d numpy array of size `q`
        Terminal solid solubility for precipitation [wt ppm].
    TSSd : 1d numpy array of size `q`
        Terminal solid solubility for dissolution [wt ppm].
    growth_fraction : float
        Fraction of precipitation corresponding to growth of hydrides [/].
    nucleation_fraction : float
        Fraction of precipitation corresponding to nucleation of hydrides [/].
    temperature_profile : 1d numpy array of size `q`
        Temperature profile in the sample [K].
    t : float
        Time [s].
    test_t : float
        Time since last output [s].
    x : 1d numpy array of size `q`
        Array with the positions determined by the mesh [m].

    """
    
    ########################## PARSE INPUT DATA ###############################
    
    # Thermal history
    
    input_temperature_positions = temperature_data[0, 1:]
    input_time_stamps = temperature_data[1:, 0]
    input_temperature = temperature_data[1:, 1:]
    
    # Experimental data
    
    slab_length, flux_left, flux_right, liner_solubility_factor, liner_width = experiment_data
    
    # Simulation data
    
    name, n_nodes, dt_div, test_t_set, criterion, model = simulation_data

    # Initial hydrogen distribution
    
    if hydrogen_data.ndim == 2:
        input_hydrogen_positions = hydrogen_data[0]
        input_hydrogen_initial_concentration = hydrogen_data[1]
        
    # Thermal history
    
    input_temperature_positions = temperature_data[0, 1:]
    input_time_stamps = temperature_data[1:, 0]
    input_temperature = temperature_data[1:, 1:]
    
    # Model parameters
    
    TSSd0 = model_parameters['TSSd0']
    Q_TSSd = model_parameters['Q_TSSd']
    TSSp0 = model_parameters['TSSp0']
    Q_TSSp = model_parameters['Q_TSSp']
    D0 = model_parameters['D0']
    Ed = model_parameters['Ed']
    Q_star = model_parameters['Q_star']
    V_star = model_parameters['V_star']
    Eg = model_parameters['Eg']
    avrami_parameter = model_parameters['avrami_parameter']
    
    Eth0 = model_parameters['Eth0']
    Eth1 = model_parameters['Eth1']
    Eth2 = model_parameters['Eth2']
    Eth3 = model_parameters['Eth3']
    
    K_G_mob0 = model_parameters['K_G_mob0']
    K_G_th0 = model_parameters['K_G_th0']
    K_N_0 = model_parameters['K_N_0']
    K_D_0 = model_parameters['K_D_0']
    
    K_G_th0 = model_parameters['K_G_th0']
    K_N_0 = model_parameters['K_N_0']
    K_D_0 = model_parameters['K_D_0']    
    
    tau = model_parameters['tau']
    g = model_parameters['g']
    delta = model_parameters['delta']
    
    ####################### CALCULATION VARIABLES #############################
    
    # Change of units to MKS (minutes to seconds, celsius to kelvin).
    
    input_time_stamps = input_time_stamps * 60
    input_temperature = input_temperature + 273
    
    # Calculate spatial discretization.
    
    dx = slab_length/n_nodes
    x = np.linspace(-dx, slab_length+dx, n_nodes+3)
    
    # Correct input positions so that each position is equal to a node position.
    
    input_temperature_positions = correct_input_positions(input_temperature_positions, x)
    
    if hydrogen_data.ndim == 2:
        input_hydrogen_positions = correct_input_positions(input_hydrogen_positions, x)
    
    # Calculation variables. The variable `test_t` is used to compare it to
    # `test_t_set` (set by the user) to know if output needs to be saved or not.
    
    t = 0
    test_t = 0
    t_max = input_time_stamps[-1]
    
    hydrostatic_stress = np.zeros(n_nodes+3)
    temperature_profile = np.zeros((n_nodes+3, 2))
    TSSp_eff = np.zeros(n_nodes+3)
    TSSd_eff = np.zeros(n_nodes+3)
    
    # Setting starting temperature and hydrogen content throghout the sample,
    # using linear interpolation for the input data. Use these as initial 
    # reference profiles.
    
    initial_temperature_profile = spatial_interpolation(x, input_temperature_positions, input_temperature[0])
    
    if hydrogen_data.ndim == 2:
        initial_hydrogen_profile = spatial_interpolation(x, input_hydrogen_positions, input_hydrogen_initial_concentration)
        ref_hydrogen_profile = initial_hydrogen_profile
    elif hydrogen_data.ndim == 1:
        ref_hydrogen_profile = np.zeros(n_nodes+3)
        
    ref_temperature_profile = initial_temperature_profile
    
    # Calculate initial TSSp and TSSd, used to find how much of the initial 
    # hydrogen is in solid solution and how much in hydrides.
    
    TSSp, TSSd = calc_terminal_solid_solubility(TSSp0, Q_TSSp, TSSd0, Q_TSSd, 
                                                initial_temperature_profile, x, 
                                                liner_width, liner_solubility_factor)
    
    # Calculate the distribution of initial hydrogen in solid solution (Css) 
    # and in hydrides (Cp), assuming equilibrium. 
    
    # If the input is a one dimensional array, it is assumed that the values
    # correspond to the initial hydrogen in solid solution and hydrogen in
    # hydrides (i. e. np.array([init_Css, init_Cp])), with a uniform hydrogen 
    # distribution (this input is used only for testing purposes).
    
    if hydrogen_data.ndim == 2:
        Css, Cp = set_initial_hydrogen_distribution(n_nodes, TSSd, initial_hydrogen_profile)
    elif hydrogen_data.ndim == 1:
        Css = hydrogen_data[0]*np.ones((n_nodes+3, 2))
        Cp = hydrogen_data[1]*np.ones(n_nodes+3)
     
    # Variables for mHNGD calculation.
    
    temperature_change_flag = np.ones(n_nodes+3, dtype='bool')
    t_0 = np.zeros(n_nodes+3)
    
    ########################## START TIME LOOP ################################
    
    while t < t_max:
        
        # Use linear interpolation to find the temperature values at the input 
        # positions in the sample for time t.
        
        temperatures_at_input_pos = time_interpolation(t, input_time_stamps, input_temperature)
        
        # Compute the temperature profile in the sample.
        
        temperature_profile[:,1] = spatial_interpolation(x, input_temperature_positions, temperatures_at_input_pos)
        
        # Beginning the precipitation/dissolution algorithm.
        
        K_N, K_G, K_D, max_K_N, max_K_G, max_K_D = calc_kinetic_parameters(Eth0, Eth1, Eth2, Eth3, 
                                                                           K_G_mob0, K_G_th0, K_D_0, K_N_0, 
                                                                           Ed, Eg, temperature_profile[:,1], 
                                                                           Cp, Css[:,0], TSSd)        
        # Recalculating the diffusion coefficient throughout the sample.
        
        diffusion_coefficient = calc_diffusion_coefficient(D0, Ed, temperature_profile[:,1])
        
        # Refreshing solubility limit values - accounting for temperature changes.
        
        TSSp, TSSd = calc_terminal_solid_solubility(TSSp0, Q_TSSp, TSSd0, Q_TSSd, 
                                                    temperature_profile[:,1], x, liner_width, 
                                                    liner_solubility_factor)
        
        if model == 'mHNGD':
            
            # Implementation of the modified HNGD model (mHNGD). Improves the
            # prediction of the thickness of the hydride rim under a temperature
            # gradient by improving the solubility limits under two main hypothesis:
            # a time-dependency to the supersolubility (TSSp) and a hydride-content
            # dependency to the solubility (TSSd). See Ref. [4] for details.
            
            diff_temp = temperature_profile[:,1] - temperature_profile[:,0]
            
            # t_0 keeps track of the time that has passed since there was no
            # temperature change (i.e. temperature has remained constant at
            # each node since t_0 = t).
            
            mask1 = np.logical_and(diff_temp == 0, temperature_change_flag)
            temperature_change_flag[mask1] = False
            t_0[mask1] = t
            
            # Effective TSSp calculation.
            
            mask2 = diff_temp == 0
            TSSp_eff[mask2] = TSSd[mask2] + (TSSp[mask2]-TSSd[mask2])*np.exp(-(t-t_0[mask2])/tau)
            TSSp_eff[~mask2] = TSSp[~mask2]
            temperature_change_flag[~mask2] = True
            t_0[~mask2] = 0
            
            # Effective TSSd calculation.
            
            Cp_at = ppm_to_atomic_fraction(Cp)
            C_delta = -9.93e-11*temperature_profile[:,1]**3 + 8.48e-8*temperature_profile[:,1]**2 - 5.73e-5*temperature_profile[:,1] + 0.623
            TSSd_at = ppm_to_atomic_fraction(TSSd)
            lever_rule_Cp = Cp_at / (C_delta - TSSd_at)
            polynomial = g*lever_rule_Cp - ((1-delta)*TSSd + g)*lever_rule_Cp**2
            TSSd_eff = TSSd + polynomial
            
            # Use effective values of TSSd and TSSp.
            
            TSSd = TSSd_eff
            TSSp = TSSp_eff
            
        # Flux diffusion equation.
        
        J_diff = calc_diffusion_flux(n_nodes, flux_left, flux_right, Q_star, V_star, 
                                     dx, temperature_profile[:,1], diffusion_coefficient, Css[:,1], 
                                     hydrostatic_stress)
        
        # Calculate reaction type in order to choose appropiate time step.
        
        reaction_type = calc_reaction_type(Css[:,1], Cp, TSSd, TSSp, activate_nucleation, activate_dissolution, activate_growth)
                    
        # Calculate time step.
        
        dt = time_step(dx, max_K_G, max_K_N, max_K_D, Q_star, dt_div, diffusion_coefficient, temperature_profile[:,1], test_t_set, reaction_type)
        
        # Solve the diffusion equation for hydrogen in solid solution.
        
        Css[1:-1, 1] = Css[1:-1, 0] - np.diff(J_diff) / dx*dt
        
        # Zero flux at the boundaries
        
        Css[0,1] = Css[1,1]
        Css[-1,1] = Css[-2,1]
        
        Cp[0] = Cp[1]
        Cp[-1] = Cp[-2]
        
        # Precipitation/dissolution equations.
          
        Css, Cp, growth, nucleation = precipitation_dissolution(n_nodes, dt, Css, Cp, TSSd, TSSp, K_N, K_G, K_D, avrami_parameter,
                                                                activate_nucleation, activate_growth, activate_dissolution)   

        # Check if temperature or hydrogen profiles changed sufficiently to 
        # warrant an output.
        
        total_hydrogen_profile = Css[:,1] + Cp
        
        hydrogen_profile_changed = profile_changed(total_hydrogen_profile, ref_hydrogen_profile, criterion)
        temperature_profile_changed = profile_changed(temperature_profile[:,1], ref_temperature_profile, criterion)
        
        if hydrogen_profile_changed:
            ref_hydrogen_profile = np.copy(total_hydrogen_profile)
        
        if temperature_profile_changed:
            ref_temperature_profile = np.copy(temperature_profile[:,1])
            
        # Saving the output at every time step could result in very large 
        # files. So, an output at this time step will be shown if: (1) either 
        # the hydrogen or temperature profiles changed significantly, or 
        # (2) the time is one of the time stamps specified in the input by the 
        # user, or (3) some time `test_t_set` has passed without an output.
       
        if hydrogen_profile_changed or temperature_profile_changed or t in input_time_stamps or test_t >= test_t_set:
            
            # Growth and nucleation fraction: which fraction of the hydride
            # precipitation corresponds to growth and which fraction to 
            # nucleation.
            
            if growth > 0 or nucleation > 0:
                growth_fraction = growth/(growth+nucleation)
                nucleation_fraction = nucleation/(growth+nucleation)
            else:
                growth_fraction = 0
                nucleation_fraction = 0
            
            # Output
            
            yield Css[:,0], Cp[:], TSSp, TSSd, growth_fraction, nucleation_fraction, temperature_profile[:,1], t, test_t, x
            
            # Resetting to start counting again
            
            test_t = 0
        
        # Actualize values for next time step
        
        Css[:,0] = Css[:,1]
        temperature_profile[:,0] = temperature_profile[:,1]
        
        # Next time step
        
        t = t + dt
        test_t = test_t + dt
