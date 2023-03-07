using Roots
using Test

"""
    topfarm_monopile_cost(depth, refTCC, nturbines)

Calculates an estimated cost of a monopile at a given depth using Equation 7 in 'TOPFARM: Multi-fidelity optimization of wind farms'
by  Réthoré et. al (2013)

# Arguments
-`depth::Float`: The water depth at the site in meters.
-`refTCC::Float`: Reference Turbine Capital Cost for the entire wind farm in USD/kW
-`nturbines::Int`: The number of turbines in the wind farm.
"""
function topfarm_monopile_cost(depth, refTCC, nturbines)  
    # CTg = 2% of the total turbine cost, CTr = 20% of the total turbine cost 
    CTr = refTCC*0.2/nturbines
    CTg = refTCC*0.02/nturbines
    return CTr + CTg*(depth-8)
end

"""
    design_monopile(mean_windspeed, site_depth, rotor_diameter, hub_height, rated_windspeed)

Calculates the required monopile design at a given site and returns the associated cost. Uses the process described in 
'10 Steps for Monopile Design' by Arany & Bhattacharya (2017) and used in NREL's  WISDEM package.

# Arguments
-`mean_windspeed::Float`: The mean wind speed in the wind farm area over a year in meters/second.
-`site_depth::Float`: The water depth at the site of the wind turbine in meters.
-`rotor_diameter::Float`: The rotor diameter of the wind turbines in meters.
-`hub_height::Float`: The hub height of the wind turbine in meters.
-`rated_windspeed::Float`: The rated windspeed of the wind turbines in meters/second.
"""
function design_monopile(mean_windspeed, site_depth, rotor_diameter, hub_height, rated_windspeed)
    yield_stress = 355000000 # Pa 
    material_factor = 1.1
    monopile_steel_cost = 700 # USD/ton gathered google search
    tp_steel_cost = 700 # general cost of steel per ton

    M_50y = calculate_50year_wind_moment(mean_windspeed, site_depth, rotor_diameter, hub_height, rated_windspeed)

    design_diameter = monopile_diameter(yield_stress, material_factor, M_50y)
    design_thickness = monopile_thickness(design_diameter)

    M_pile = calculate_pile_moment(design_diameter, design_thickness)

    air_gap = 10 #m
    design_length = air_gap + site_depth + pile_embedment_length(M_pile) 

    mass = monopile_mass(design_diameter, design_thickness, design_length)

    design_cost = mass * monopile_steel_cost

    # Include cost for transition piece based on monopile design
    connection_thickness = 0.0 # bolted connection
    L_tp = 25 # m
    dens_tp = 7860 # kg/m^3
    mass_tp = (dens_tp * (design_diameter + 2 * connection_thickness + design_thickness) * pi * design_thickness * L_tp) / 907.185 # convert to tons
    tp_cost = mass_tp * tp_steel_cost

    return (design_cost + tp_cost)
end

"""
    monopile_diameter(yield_stress, material_factor, M_50y)

Calculates the required monopile diameter. Uses Equations 99 and 101 from from '10 Steps for Monopile Design' by Arany & Bhattacharya (2017).

# Arguments
-`yield_stress::Float`: The yield stress of the monopile in Pascals.
-`material_factor::Float`: The material factor.
-`M_50y::Float`: The maximum moment on the monopile over a 50 year period in Newton-meters.
"""
function monopile_diameter(yield_stress, material_factor, M_50y)
    A = (yield_stress * pi) / (4 * material_factor * M_50y)
    f(Dp) = A * ((0.99 * Dp - 0.00635) ^ 3) * (0.00635 + 0.01 * Dp) - Dp

    # Find the root of equation f to identify a valid diameter.
    design_diameter = find_zero(f, 10)

    return design_diameter
end

"""
    calculate_50year_wind_moment(mean_windspeed, site_depth, rotor_diameter, hub_height,rated_windspeed)

Calculates the maximum moment on the monopile from the wind during a 50 year period. Uses Equation 30
from '10 Steps for Monopile Design' by Arany & Bhattacharya (2017).

# Arguments
-`mean_windspeed::Float`: The mean wind speed in the wind farm area over a year in meters/second.
-`site_depth::Float`: The water depth at the site of the wind turbine in meters.
-`rotor_diameter::Float`: The rotor diameter of the wind turbines in meters.
-`hub_height::Float`: The hub height of the wind turbine in meters.
-`rated_windspeed::Float`: The rated windspeed of the wind turbines in meters/second.
"""
function calculate_50year_wind_moment(mean_windspeed, site_depth, rotor_diameter, hub_height, rated_windspeed)
    load_factor = 3.375
    F_50y = calculate_50year_wind_load(mean_windspeed, rotor_diameter, rated_windspeed)

    M_50y = F_50y * (site_depth + hub_height)

    return M_50y * load_factor
end

"""
    calculate_50year_wind_load(mean_windspeed, rotor_diameter, rated_windspeed)

Calculates the maximum load on the monopile due to the wind over 50 years. Uses Equation 29
from '10 Steps for Monopile Design' by Arany & Bhattacharya (2017).

# Arguments
-`mean_windspeed::Float`: The mean wind speed in the wind farm area over a year in meters/second.
-`rotor_diameter::Float`: The rotor diameter of the wind turbines in meters.
-`rated_windspeed::Float`: The rated windspeed of the wind turbines in meters/second.
"""
function calculate_50year_wind_load(mean_windspeed, rotor_diameter, rated_windspeed)
    air_density = 1.225
    swept_area = pi * (rotor_diameter / 2)^2

    coeff_thrust = calculate_thrust_coefficient(rated_windspeed)

    U_extreme = calculate_50year_extreme_gust(mean_windspeed, rotor_diameter, rated_windspeed)

    F_50y = 0.5 * air_density * swept_area * coeff_thrust * ((rated_windspeed + U_extreme)^2)

    return F_50y
end

"""
    calculate_thrust_coefficient(rated_windspeed)

Calculates the thrust coefficient of the wind turbines.

# Arguments
-`rated_windspeed::Float`: The rated wind speed of the turbines in the wind farm in meters/second.
"""
function calculate_thrust_coefficient(rated_windspeed)
    ct = min(3.5 * (2 * rated_windspeed + 3.5) / (rated_windspeed^2), 1)

    return ct
end

"""
    pile_embedment_length(Mp)

Calculates the embedment length of the monopile using equation 7 from '10 Steps for Monopile Design'
by Arany & Bhattacharya (2017). Uses assumed values for the monopile's modulus of elasticity and 
soil coefficient.

# Arguments:
-`Mp::Float`: Moment of the monopile in Newton-meters.
"""
function pile_embedment_length(Mp)
    # TODO: Add inputs for main COE file and variables in params to specify monopile
    # modulus, soil coefficient, and other variables (cost of steel, etc)
    monopile_modulus = 200e9 # Pa
    soil_coefficient = 4000000 # N/m^3

    Length = 2 * ((monopile_modulus * Mp) / soil_coefficient)^0.2

    return Length
end

"""
    calculate_50year_extreme_gust(mean_windspeed, rotor_diameter, rated_windspeed)

Calculates the estimated extreme wind gust over 50 years. Uses Equation 28
from '10 Steps for Monopile Design' by Arany & Bhattacharya (2017).

# Abstract
-`mean_windspeed::Float`: The mean wind speed in the wind farm area over a year in meters/second.
-`rotor_diameter::Float`: The rotor diameter of the wind turbines in meters.
-`rated_windspeed::Float`: The rated windspeed of the wind turbines in meters/second.
"""
function calculate_50year_extreme_gust(mean_windspeed, rotor_diameter, rated_windspeed)
    length_scale = 340.2

    U_50y = calculate_50year_extreme_ws(mean_windspeed)
    U_1y = 0.8 * U_50y

    U_eog = min(
            (1.35 * (U_1y - rated_windspeed)),
            (3.3 * 0.11 * U_1y) / (1 + (0.1 * rotor_diameter) / (length_scale / 8))
    )

    return U_eog
end

"""
    calculate_50year_extreme_ws(mean_windspeed)

Calculates the estimated extreme wind speed over 50 years. Uses Equation 27
from '10 Steps for Monopile Design' by Arany & Bhattacharya (2017).

# Abstract
-`mean_windspeed::Float`: The mean wind speed in the wind farm area over a year in meters/second.
"""
function calculate_50year_extreme_ws(mean_windspeed)
    scale_factor = mean_windspeed
    shape_factor = 2
    U_50y = scale_factor * (-log(1 - 0.98 ^ (1 / 52596))) ^ (1 / shape_factor)

    return U_50y
end

"""
    monopile_mass(D, t, L)

Calculates the mass of a monopile in US tons. Assumes the density of the monopile steel is 7,860 kg/m^3

# Arguments
-`D::Float`: Diameter of the monopile in meters
-`t::Float`: Thickness of the monopile in meters
-`L::Float`: Length of the monopile in meters
"""
function monopile_mass(D, t, L)
    # Line 328 of monopile_design.py in WISDEM from NREL
    density = 7860 # kg/m3
    volume = (pi / 4) * (D^2 - (D - t)^2) * L
    mass = density * volume / 907.185 # convert the mass to tons

    return mass # In tons
end

"""
    monopile_thickness(D)

Calculates the thickness of a monopile using Equation 1 from '10 Steps for Monopile Design'
by Arany & Bhattacharya (2017)

# Arguments
-`D::Float`: Diameter of monopile in meters
"""
function monopile_thickness(D)
    thickness = 0.00635 + D/100

    return thickness
end

"""
    calculate_pile_moment(D, t)

Calculates the moment on a monopile of a given diameter and thickness

# Arguments
-`D::Float`: Diameter of monopile in meters
-`t::Float`: Thickness of monopile in meters
"""
function calculate_pile_moment(D, t)
    Ip = 0.125 * ((D - t)^3) * t * pi

    return Ip
end