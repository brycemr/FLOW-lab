using FLOWFarm; const ff = FLOWFarm
using PyPlot; const plt = PyPlot
using VectorizedRoutines.Matlab: meshgrid
using SNOW
using Plots
using Distributions

include("../FLOW lab/Depths.jl")
include("../FLOW lab/Monopile.jl")

# set wind farm boundary parameters in meters
boundarycenter = [0.0, 0.0]
boundaryradius = hypot(707,707)

#--------- SET DEPTH SETTINGS -------------#

gaus1D = false
gaus2D = true
gausDepthOutward = true
linDepth = false
gausMinDepth = 10
gausMaxDepth = 60
linStartDepth = 23.5
linDepthVariance = 0
monopileCostModel = "Mass" # Select between Linear and Mass

depthBoundaries = [-boundaryradius*1.1 boundaryradius*1.1; -boundaryradius*1.1 boundaryradius*1.1]

if gaus2D
    gausStartDepth = gausMinDepth
    means = [0.0, 0.0] 
    coVariants = [9.0e5 6.0e5; 6.0e5 9.0e5] 
    distribution = MvNormal(means, coVariants)
else
    gausStartDepth = gausMinDepth
    depthStdDev = (gausMaxDepth - gausMinDepth)/3
    distribution = Normal(gausStartDepth, depthStdDev)
end

turbinex = [-500.0, -500.0, -500.0, 0.0, 0.0, 0.0, 500.0, 500.0, 500.0]
turbiney = [-500.0, 0.0, 500.0, -500.0, 0.0, 500.0, -500.0, 0.0, 500.0]

nturbines = length(turbinex)

turbinez = zeros(nturbines)

turbineyaw = zeros(nturbines)

# set wind turbine design parameters
rotordiameter = zeros(nturbines) .+ 125 # m 
hubheight = zeros(nturbines) .+ 90      # m
cutinspeed = zeros(nturbines) .+ 3.0    # m/s
cutoutspeed = zeros(nturbines) .+ 25.0   # m/s 
ratedspeed = zeros(nturbines) .+ 11.4   # m/s 
ratedpower = zeros(nturbines) .+ 5.0E6  # W (5 MW)
generatorefficiency = ones(nturbines)

#---------- VISUALIZING THE WIND FARM LAYOUT ----------------#

fig, ax1 = plt.subplots(1)

ff.plotlayout!(ax1, turbinex, turbiney, rotordiameter)

ax1.set(xlabel="Easting (m)", ylabel="Northing (m)")

circle = matplotlib.patches.Circle((0.0, 0.0), boundaryradius, fill=false, color="k")
ax1.add_patch(circle)

ax1.set(xlim=[-boundaryradius, boundaryradius].*1.01, ylim=[-boundaryradius, boundaryradius].*1.01)

mapDepths(50, [-boundaryradius*1.1, boundaryradius*1.1], [-boundaryradius*1.1, boundaryradius*1.1], distribution, gausStartDepth, linStartDepth, linDepthVariance, gausDepthOutward, linDepth, gaus1D, gaus2D, gausMaxDepth, gausMinDepth)

display(fig)

# get sample points for rotor swept area for determining inflow wind speed
nsamplepoints = 50
rotorsamplepointsy, rotorsamplepointsz = ff.rotor_sample_points(nsamplepoints, method="sunflower")

#--------- SETUP WIND RESOURCE ----------#
windspeed = 9.74 # m/s 
airdensity = 1.1716 # kg/m^3
ambientti = 0.1 # %
shearexponent = 0.12
ndirections = 5
winddirections = collect(range(0, 2*pi*(1-1/ndirections), length=ndirections))
windspeeds = ones(ndirections).*windspeed
windprobabilities = ones(ndirections).*(1.0/ndirections)
ambienttis = ones(ndirections).*ambientti
measurementheight = ones(ndirections).*hubheight[1]

# initialize the wind shear model
windshearmodel = ff.PowerLawWindShear(shearexponent)

# initialize the wind resource definition
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities,
measurementheight, airdensity, ambienttis, windshearmodel)

# visualize the wind resource
fig = ff.plotwindresource!(windresource)
display(fig)
#--------- SET INITIAL TURBINE COSTS------------#

# Values gathered from NREL 2020 Cost of Energy report, BOS has foundation cost removed
#refTCC = 1301.0 # $/kW
refTCC = 1952.0
#refBOS = 1782.0 - 474.0 # $/kW
refBOS = 4420 - 474.0
refFC = 672.0 # # $/kW
refFCR = 0.0764 # /year
#refCapEx = 110.0 # $/kw/year
refCapEx = 140.56
refLCOE = ff.Levelized(refTCC, refBOS, refFC, refFCR, refCapEx)

#--------- SETTING UP ANALYSIS MODELS ----------#
powermodel = ff.PowerModelPowerCurveCubic()

powermodels = Vector{typeof(powermodel)}(undef, nturbines)
for i = 1:nturbines
    powermodels[i] = powermodel
end

ctmodel = ff.ThrustModelConstantCt(0.65)
ctmodels = Vector{typeof(ctmodel)}(undef, nturbines)
for i = 1:nturbines
    ctmodels[i] = ctmodel
end

wakedeficitmodel = ff.JensenCosine()
wakedeflectionmodel = ff.GaussYawDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelMaxTI()
modelset = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

# --------SETTING UP CONSTRAINTS AND OPTIMIZATION FUNCTIONS------------#

# scale objective derivatives to be between 0 and 1
objectivescale = 1.0E-3

# scale boundary constraint derivatives to be between 0 and 1
constraintscaleboundary = 1.0E-5

# scale spacing constraint derivatives to be between 0 and 1
constraintscalespacing = 1.0

# set the minimum spacing between turbines 
minimumspacing = 159.0*1.2

# set up a struct for use in optimization functions, these are the non-differentiated parameters
mutable struct params_struct{}
    modelset
    rotorsamplepointsy
    rotorsamplepointsz
    turbinez
    ambientti
    rotordiameter
    boundarycenter
    boundaryradius
    objectivescale
    constraintscaleboundary
    constraintscalespacing
    minimumspacing
    hubheight
    turbineyaw
    ctmodels
    generatorefficiency
    cutinspeed
    cutoutspeed
    ratedspeed
    ratedpower
    windresource
    powermodels
    gaus1D
    gaus2D
    gausDepthOutward
    linDepth
    linDepthVariance
    gausMinDepth
    gausMaxDepth
    linStartDepth
    distribution
    depthBoundaries
    refLCOE
    meanWindspeed
    monopileCostModel
end

params = params_struct(modelset, rotorsamplepointsy, rotorsamplepointsz, turbinez, ambientti,
    rotordiameter, boundarycenter, boundaryradius, objectivescale, constraintscaleboundary,
    constraintscalespacing, minimumspacing, hubheight, turbineyaw,
    ctmodels, generatorefficiency, cutinspeed, cutoutspeed, ratedspeed, ratedpower,
    windresource, powermodels, gaus1D, gaus2D, gausDepthOutward, linDepth, linDepthVariance,
     gausMinDepth, gausMaxDepth, linStartDepth, distribution, depthBoundaries, refLCOE, windspeed, monopileCostModel)


# Set up wrapper functions for the objective and constraints
function boundary_wrapper(x, params)
    # include the relevant params
    boundarycenter = params.boundarycenter
    boundaryradius = params.boundaryradius
    constraintscaleboundary = params.constraintscaleboundary

    nturbines = Int(length(x)/2)

    # extract the x and y locations of turbines from design variables Vector
    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]

    # get and return boundary distances
    return ff.circle_boundary(boundarycenter, boundaryradius, turbinex, turbiney).*constraintscaleboundary
end

# Set up spacing constraint wrapper function
function spacing_wrapper(x, params)
    # include relevant params
    rotordiameter = params.rotordiameter
    constraintscalespacing = params.constraintscalespacing
    minimumspacing = params.minimumspacing

    nturbines = Int(length(x)/2)

    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]

    return constraintscalespacing.*(minimumspacing .- ff.turbine_spacing(turbinex, turbiney))
end

function total_monopile_cost_wrapper(x, params)
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
        
    end

    return cost
end

# Set up the aep wrapper function
function coe_wrapper(x, params)
    turbinez = params.turbinez
    rotordiameter = params.rotordiameter
    hubheight = params.hubheight
    turbineyaw =params.turbineyaw
    ctmodels = params.ctmodels
    generatorefficiency = params.generatorefficiency
    cutinspeed = params.cutinspeed
    cutoutspeed = params.cutoutspeed
    ratedspeed = params.ratedspeed
    ratedpower = params.ratedpower
    windresource = params.windresource
    powermodels = params.powermodels
    modelset = params.modelset
    rotorsamplepointsy = params.rotorsamplepointsy
    rotorsamplepointsz = params.rotorsamplepointsy
    objectivescale = params.objectivescale
    refCost = params.refLCOE

    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]

    # calculate_aep
    aep = ff.calculate_aep(turbinex, turbiney, turbinez, rotordiameter,
    hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
    cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
    rotor_sample_points_y=rotorsamplepointsy,rotor_sample_points_z=rotorsamplepointsz)

    monopilesCost = total_monopile_cost(x, params)

    cost = ff.Levelized(refCost.TCC, refCost.BOS+monopilesCost, refCost.FC,
    refCost.FCR, refCost.OpEx)    
    
    total_ratedpower = sum(ratedpower)
    coe_aep = aep/total_ratedpower # For the LCOE equation we need aep in MWh/MW/year not Wh/year (same as hr/year)
    coe_ratedpower = ratedpower./1000 # needs to be in units of kw
    coe = objectivescale.*ff.cost_of_energy(rotordiameter, hubheight, coe_ratedpower,  coe_aep, cost)

    return coe
end

# Set up the optimization problem wrapper function
function wind_farm_opt!(g, x, params)
    nturbines = Int(length(x)/2)

    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x, params)

    # calculate boundary constraint and jacobian
    boundary_con = boundary_wrapper(x, params)

    # combine constraint values and jacobians into overall constraint value and jacobian
    g[1:(end-nturbines)] = spacing_con[:]
    g[end-nturbines+1:end] = boundary_con[:]

    obj = coe_wrapper(x, params)[1]
    #print("\n OBJECTIVE \n")
    #println(obj)
    return obj
end

# -------------------PRINT INITIAL VALUES-------------------------------
# Annual Energy production aep
aep = ff.calculate_aep(turbinex, turbiney, turbinez, rotordiameter,
    hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
    cutoutspeed, ratedspeed, ratedpower, windresource, powermodels,
    modelset, rotor_sample_points_y=rotorsamplepointsy, 
    rotor_sample_points_z=rotorsamplepointsz)

# Cost of Energy
currPos = [copy(turbinex);copy(turbiney)]
monopilesCost = total_monopile_cost(currPos, params)
cost = ff.Levelized(refLCOE.TCC+monopilesCost, refLCOE.BOS, refLCOE.FC,
    refLCOE.FCR, refLCOE.OpEx)    

total_ratedpower = sum(ratedpower)
coe_aep = aep/total_ratedpower # For the LCOE equation we need aep in MWh/MW/year not Wh/year (same as hr/year)
coe_ratedpower = ratedpower./1000 # needs to be in units of kw
coe_initial = ff.cost_of_energy(rotordiameter, hubheight, coe_ratedpower, coe_aep, cost)

println("$aep Watt-hours per year")
# AEP in each direction
state_aeps = ff.calculate_state_aeps(turbinex, turbiney, turbinez, rotordiameter,
        hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
        cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
        rotor_sample_points_y=rotorsamplepointsy, rotor_sample_points_z=rotorsamplepointsz,
        hours_per_year=365.25*24.0, weighted=true)

println("$state_aeps Watt-hours per year")

# generate objective function wrapper
obj_func!(g, x) = wind_farm_opt!(g, x, params)


# ------------SETUP THE OPTIMIZER --------------------------------------

# initialize design variable vector
x0 = [copy(turbinex);copy(turbiney)]

# set general lower and upper bounds for design variable
lx = zeros(length(x0)) .- boundaryradius
ux = zeros(length(x0)) .+ boundaryradius

# set general lower and upper bounds for constraints
ng = Int(nturbines + (nturbines)*(nturbines-1)/2)
lg = [-Inf*ones(Int((nturbines)*(nturbines - 1)/2)); -Inf*ones(nturbines)]
ug = [zeros(Int((nturbines)*(nturbines - 1)/2)); zeros(nturbines)]

# IPOPT options
ip_options = Dict(
    "max_iter" => 50,
    "tol" => 1e-6
)
solver = IPOPT(ip_options)

# initialize SNOW options
options = Options(solver=solver, derivatives=ForwardAD())

# -------- RUN OPTIMIZER -----------------------------
t1 = time()
xopt, fopt, info, out = minimize(obj_func!, x0, ng, lx, ux, lg, ug, options)
t2 = time()
clk = t2-t1

coefinal = fopt/objectivescale

println("Finished in : ", clk, " (s)")
println("info: ", info)
print("Init LCOE WITHOUT Monopiles: ")
println(ff.cost_of_energy(rotordiameter, hubheight, coe_ratedpower, coe_aep, refLCOE), " \$/MWh")
println("Init LCOE of Monopiles: ", monopilesCost*refLCOE.FCR/(coe_aep/1000), " \$/MWh")
println("Final LCOE of Monopiles: ", total_monopile_cost(xopt, params)*refLCOE.FCR/(coe_aep/1000), " \$/Mwh" )
println("Initial COE: ", coe_initial, " \$/MWh")
println("Final COE: ", coefinal, " \$/MWh")
println("COE improvement (%) = ", -100*(coefinal - coe_initial)/coe_initial)

# final turbine locations
turbinexopt = copy(xopt[1:nturbines])
turbineyopt = copy(xopt[nturbines+1:end])


#-----PLOT OPTIMIZED LAYOUT------------#
fig, ax = plt.subplots(1)

ff.plotlayout!(ax, turbinexopt, turbineyopt, rotordiameter)

ax.set(xlabel="Easting (m)", ylabel="Northing (m)")

# and the wind farm boundary
circle = matplotlib.patches.Circle((0.0, 0.0), boundaryradius, fill=false, color="k")
ax.add_patch(circle)

# set limits on the plot region
ax.set(xlim=[-boundaryradius, boundaryradius].*1.01, ylim=[-boundaryradius, boundaryradius].*1.01)

display(fig)