
using Roots
using Test

function total_monopile_cost(x, params)
    cost = 0
    
    nturbines = Int(length(x)/2)
    
    mean_windspeed = params.meanWindspeed
    rotordiameter = params.rotordiameter
    hub_height = params.hubheight
    rated_windspeed = params.ratedspeed
    rated_power = params.ratedpower
    depthBoundaries = params.depthBoundaries
    gaus1D = params.gaus1D
    gaus2D = params.gaus2D
    centerIsShallow = params.gausDepthOutward
    linDepth = params.linDepth
    linDepthVariance = params.linDepthVariance
    gausMaxDepth = params.gausMaxDepth
    gausMinDepth = params.gausMinDepth
    linStartDepth = params.linStartDepth
    refCost = params.refLCOE
    data = params.distribution
    monopileModel = params.monopileCostModel

    # extract x and y locations of turbines from design variables vector
    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]

    xRange = depthBoundaries[1,2] - depthBoundaries[1,1]
    yRange = depthBoundaries[2,2] - depthBoundaries[2,1]

    cost = 0

    for i = eachindex(turbinex)
        if gaus2D
            pos = [Int(ceil(1000*turbinex[i]))/1000, Int(ceil(1000*turbiney[i]))/1000]
            depth = gaussianTwoD(distribution, gausMaxDepth, gausMinDepth, pos, centerIsShallow)
        else
            linearQuant = Int(ceil(1000*(turbiney[i] - depthBoundaries[2,1] + 0.0001)/yRange))/1000
            gaussQuant = Int(ceil(1000*(turbinex[i] - depthBoundaries[1,1] + 0.0001)/xRange))/1000
            distDepth = quantile(data, gaussQuant)[1]
            difference = abs(gausMaxDepth - distDepth)
            
            depth = (gaus1D ? gaussianOneD(gausStartDepth, difference, gausDepthOutward) : 0)
            depth = (linDepth ? depth + linearDepth(linDepthVariance, linStartDepth, linearQuant) : depth)
        end
        if monopileModel == "Linear"
            latest = topfarm_monopile_cost(depth, refCost.TCC, nturbines)
            cost += latest
        else
            latest = design_monopile(mean_windspeed, depth, rotordiameter[i], hub_height[i], rated_windspeed[i])
            cost += latest / (rated_power[i]/1000) # Convert from USD to USD/kW for use in LCOE equation
        end
        #println("Latest:$latest Cost:$cost Depth:$depth mean_windspeed:$mean_windspeed")
    end

    return cost
end

# This function calculates the cost of a single monopile in a windfarm
function topfarm_monopile_cost(depth, refTCC, nturbines)
    # Using equation 7 from TOPFARM: Multi-fidelity optimization    
    # CTg = 2% of the total turbine cost, CTr = 20% of the total turbine cost 
    CTr = refTCC*0.2/nturbines
    CTg = refTCC*0.02/nturbines
    return CTr + CTg*(depth-8)
end

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

function monopile_diameter(yield_stress, material_factor, M_50y)
    # Equations for determining diameter, Eq 99 and 101
    A = (yield_stress * pi) / (4 * material_factor * M_50y)
    f(Dp) = A * ((0.99 * Dp - 0.00635) ^ 3) * (0.00635 + 0.01 * Dp) - Dp

    design_diameter = find_zero(f, 10)

    return design_diameter
end

function calculate_50year_wind_moment(mean_windspeed, site_depth, rotor_diameter, hub_height, rated_windspeed)
    # Equation 30 from Arany & Bhattacharya (2017)
    load_factor = 3.375
    F_50y = calculate_50year_wind_load(mean_windspeed, rotor_diameter, rated_windspeed)

    M_50y = F_50y * (site_depth + hub_height)

    return M_50y * load_factor
end

function calculate_50year_wind_load(mean_windspeed, rotor_diameter, rated_windspeed)
    # Equation 29 from Arany & Bhattacharya (2017)
    air_density = 1.225
    swept_area = pi * (rotor_diameter / 2)^2

    coeff_thrust = calculate_thrust_coefficient(rated_windspeed)

    U_extreme = calculate_50year_extreme_gust(mean_windspeed, rotor_diameter, rated_windspeed)

    F_50y = 0.5 * air_density * swept_area * coeff_thrust * ((rated_windspeed + U_extreme)^2)

    return F_50y
end

function calculate_thrust_coefficient(rated_windspeed)
    ct = min(3.5 * (2 * rated_windspeed + 3.5) / (rated_windspeed^2), 1)

    return ct
end

function pile_embedment_length(Mp)
    # Uses eq 7
    # TODO: Add inputs for main COE file and variables in params to specify monopile
    # modulus, soil coefficient, and other variables (cost of steel, etc)
    monopile_modulus = 200e9 # Pa
    soil_coefficient = 4000000 # N/m^3

    Length = 2 * ((monopile_modulus * Mp) / soil_coefficient)^0.2

    return Length
end

function calculate_50year_extreme_gust(mean_windspeed, rotor_diameter, rated_windspeed)
    # Equation 28 

    length_scale = 340.2

    U_50y = calculate_50year_extreme_ws(mean_windspeed)
    U_1y = 0.8 * U_50y

    U_eog = min(
            (1.35 * (U_1y - rated_windspeed)),
            (3.3 * 0.11 * U_1y) / (1 + (0.1 * rotor_diameter) / (length_scale / 8))
    )

    return U_eog
end

function calculate_50year_extreme_ws(mean_windspeed)
    # Equation 27
    scale_factor = mean_windspeed
    shape_factor = 2
    U_50y = scale_factor * (-log(1 - 0.98 ^ (1 / 52596))) ^ (1 / shape_factor)

    return U_50y
end

function monopile_mass(D, t, L)
    # Line 328 of monopile_design.py in WISDEM from NREL
    density = 7860 # kg/m3
    volume = (pi / 4) * (D^2 - (D - t)^2) * L
    mass = density * volume / 907.185 # Why 907.185

    return mass
end

function monopile_thickness(D)
    # Equation 1 from Arany & Bhattacharya (2017) 10 steps for monopile design
    thickness = 0.00635 + D/100

    return thickness
end

function calculate_pile_moment(D, t)
    Ip = 0.125 * ((D - t)^3) * t * pi

    return Ip
end

@testset "50 year extreme Wind Speed" begin
    @test calculate_50year_extreme_ws(9.74) == 37.43548973876968
    @test calculate_50year_extreme_ws(7) == 26.904356075091144
    @test calculate_50year_extreme_ws(8) == 30.74783551438988
end

@testset "50 year extreme Wind Gust" begin
    @test calculate_50year_extreme_gust(9.74, 125, 11.4) == 8.401646451820062
    @test calculate_50year_extreme_gust(7, 125, 11.4) == 6.038144267221809
    @test calculate_50year_extreme_gust(7, 100, 11.4) == 6.325538092410854
    @test calculate_50year_extreme_gust(7, 125, 10) == 6.038144267221809
end

@testset "Thrust Coefficient" begin
    @test calculate_thrust_coefficient(11.4) == 0.7082948599569098
    @test calculate_thrust_coefficient(7) == 1
    @test calculate_thrust_coefficient(8) == 1
    @test calculate_thrust_coefficient(10) == 0.8225
end

@testset "50 year Wind Load" begin
    @test calculate_50year_wind_load(9.74, 125, 11.4) == 2087529.8529107259
    @test calculate_50year_wind_load(7, 125, 11.4) == 1618939.5140525973
    @test calculate_50year_wind_load(7, 100, 11.4) == 1070554.8426587582
    @test calculate_50year_wind_load(7, 125, 10) == 1590230.7187345193
end

@testset "50 year Wind Moment" begin
    @test calculate_50year_wind_moment(9.74, 23.5, 125, 90, 11.4) == 799654404.2806149
    @test calculate_50year_wind_moment(7, 23.5, 125, 90, 11.4) == 620155017.6017731
    @test calculate_50year_wind_moment(8, 26, 125, 90, 10) == 691345156.9397321
    @test calculate_50year_wind_moment(8, 30, 125, 80, 11.4) == 661962945.240009
end

@testset "Monopile Diameter" begin
    @test isapprox(monopile_diameter(355000000, 1.1, 799654404.2806149), 6.677656219964458; atol = 0.001)
    @test isapprox(monopile_diameter(355000000, 1.1, 620155017.6017731), 6.119550254228244; atol = 0.001)
    @test isapprox(monopile_diameter(355000000, 1.1, 691345156.9397321), 6.352343360869822; atol = 0.001)
    @test isapprox(monopile_diameter(355000000, 1.1, 661962945.240009), 6.25829318500689; atol = 0.001)
end

@testset "Monopile Thickness" begin
    @test isapprox(monopile_thickness(6.677656219964458), 0.07312656219964457; 0.001)
    @test isapprox(monopile_thickness(5), 0.056350000000000004; 0.001)
    @test isapprox(monopile_thickness(7), 0.07635; 0.001)
    @test isapprox(monopile_thickness(6.119550254228244), 0.06754550254228243; 0.001)
end

@testset "Monopile Moment" begin
    @test isapprox(calculate_pile_moment(6.677656219964458, 0.07312656219964457), 8.27295623552583; 0.001)
    @test isapprox(calculate_pile_moment(5, 0.056350000000000004), 2.6736032113211996; 0.001)
    @test isapprox(calculate_pile_moment(7, 0.07635), 9.95117225211555; 0.001)
    @test isapprox(calculate_pile_moment(6.119550254228244, 0.06754550254228243), 5.879685598856954; 0.001)
end

@testset "Embedment Length" begin
    @test isapprox(pile_embedment_length(8.27295623552583), 26.56783356260234; 0.001)
    @test isapprox(pile_embedment_length(2.6736032113211996), 21.195486401869456; 0.001)
    @test isapprox(pile_embedment_length(9.95117225211555), 27.567592805981995; 0.001)
    @test isapprox(pile_embedment_length(5.879685598856954), 24.813887979958285; 0.001)
end

@testset "Monopile Mass" begin
    @test isapprox(monopile_mass(6.677656219964458, 0.07312656219964457, (23.5+10+26.56783356260234)), 397.01164557176816; 0.001)
    @test isapprox(monopile_mass(6.119550254228244, 0.06754550254228243, (23.5+10+24.813887979958285)), 326.2353875559374; 0.001)
end
