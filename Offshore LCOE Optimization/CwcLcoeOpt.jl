using FLOWFarm; const ff = FLOWFarm
using PyPlot; const plt = PyPlot
using VectorizedRoutines.Matlab: meshgrid
using SNOW
using Plots
using Distributions

include("Depths.jl")
include("Monopile.jl")


#--------- SET DEPTH SETTINGS -------------#
boundary_min_max_latitude = [28.8299, 28.9111]
depth_map_boundary_latitude = [28.8299*0.9999, 28.9111*1.0001]

boundary_min_max_longitude = [-90.3347, -90.2414]
depth_map_boundary_longitude = [-90.3347*1.0001, -90.2414*0.9999]

boundary_center_cor = [-90.2287, 28.8299]

boundary_center_x_m = longitude_to_meters(boundary_center_cor[2], boundary_center_cor[1] - boundary_min_max_longitude[1])
boundary_center_y_m = latitude_to_meters(boundary_center_cor[1] - boundary_min_max_latitude[1])

xRange = longitude_to_meters((boundary_min_max_latitude[2] + boundary_min_max_latitude[1])/2, boundary_min_max_longitude[2] - boundary_min_max_longitude[1])
yRange = latitude_to_meters(boundary_min_max_latitude[2] - boundary_min_max_latitude[1])

depthXRange = longitude_to_meters((depth_map_boundary_latitude[2] + depth_map_boundary_latitude[1])/2, depth_map_boundary_longitude[2] - depth_map_boundary_longitude[1])
depthYRange = latitude_to_meters(depth_map_boundary_latitude[2] - depth_map_boundary_latitude[1])
depthBoundaries = [0 depthXRange; 0 depthYRange]

generateDepthPattern = false
gaus1D = true
gaus2D = false
gausDepthOutward = true
linDepth = false
mean_depth = 25
variance = mean_depth*1.0

gausMinDepth = mean_depth
gausMaxDepth = mean_depth + variance

linDepthVariance = variance
linStartDepth = mean_depth - linDepthVariance/2

if gaus2D
    gausStartDepth = gausMinDepth
    means = [0.0, 0.0]
    xxVar = 5.0e5
    xyVar = -3.5e5 
    coVariants = [xxVar xyVar; xyVar xxVar] 
    distribution = MvNormal(means, coVariants)
else
    gausStartDepth = gausMinDepth
    depthStdDev = (gausMaxDepth - gausMinDepth)/3
    distribution = Normal(gausStartDepth, depthStdDev)
end

resolution_x = 50
resolution_y = 50
depth_file_path = "C:\\Users\\bryce\\Desktop\\FLOW-lab\\Offshore LCOE Optimization\\H06155.xyz"


monopileCostModel = "Mass" # Select between Linear and Mass

# set up wind farm boundary parameters in meters
boundarycenter = (xRange/2, yRange/2)
boundaryradius = xRange
boundary_vertices = [0 0; 0 yRange; xRange yRange; xRange boundary_center_y_m; boundary_center_x_m boundary_center_y_m; boundary_center_x_m 0]
boundary_normals = [-1 0; 0 1; 1 0; 0 -1; 1 0; 0 -1]

#------------ TURBINE SETTINGS ----------------#
num_turbines = 10
turbinex = zeros(num_turbines)
turbiney = zeros(num_turbines)
turbinex .+= rand.((0:xRange/2, ))
turbiney .+= rand.((0:yRange/2, ))

# Data for best layout found after 1000 iterations at constant depth
nturbines = length(turbinex)

turbinez = zeros(nturbines)

turbineyaw = zeros(nturbines)

# set wind turbine design parameters
rotordiameter = zeros(nturbines) .+ 220 # m 
hubheight = zeros(nturbines) .+ 260      # m
cutinspeed = zeros(nturbines) .+ 3.0    # m/s
cutoutspeed = zeros(nturbines) .+ 25.0   # m/s 
ratedspeed = zeros(nturbines) .+ 11.4   # m/s 
ratedpower = zeros(nturbines) .+ 12.0E6  # W (5 MW)
generatorefficiency = ones(nturbines)

#---------- VISUALIZING THE WIND FARM LAYOUT ----------------#

fig, ax1 = plt.subplots(1)

ff.plotlayout!(ax1, turbinex, turbiney, rotordiameter)

ax1.set(xlabel="Easting (m)", ylabel="Northing (m)")

#circle = matplotlib.patches.Circle((0.0, 0.0), boundaryradius, fill=false, color="k")
square = matplotlib.patches.Rectangle((0, 0), xRange, yRange, fill=false, color="k")
ax1.add_patch(square)

ax1.set(xlim=[-xRange*0.01, xRange].*1.01, ylim=[-yRange*0.01, yRange].*1.01)

if generateDepthPattern
    mapDepths(50, [-boundaryradius*1.01, boundaryradius*1.01], [-boundaryradius*1.01, boundaryradius*1.01], distribution, gausStartDepth, linStartDepth, linDepthVariance, gausDepthOutward, linDepth, gaus1D, gaus2D, gausMaxDepth, gausMinDepth)
else
    depthMap = depthMap = import_depth_from_xyz(depth_file_path, resolution_x, resolution_y, depth_map_boundary_latitude, depth_map_boundary_longitude)
end

display(fig)

# get sample points for rotor swept area for determining inflow wind speed
nsamplepoints = 50
rotorsamplepointsy, rotorsamplepointsz = ff.rotor_sample_points(nsamplepoints, method="sunflower")

#--------- SETUP WIND RESOURCE ----------#
windspeed = 9.74 # m/s 
airdensity = 1.1716 # kg/m^3
ambientti = 0.1 # %
shearexponent = 0.12
ndirections = 12
winddirections = collect(range(0, 2*pi*(1-1/ndirections), length=ndirections))
windspeeds = ones(ndirections).*windspeed
#windprobabilities = ones(ndirections).*(1.0/ndirections)
windprobabilities = [8, 9, 7, 3, 2, 5, 9, 12, 16, 12, 10, 7] ./100
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
mutable struct params_struct5{}
    modelset
    rotorsamplepointsy
    rotorsamplepointsz
    turbinez
    ambientti
    rotordiameter
    boundaryradius
    boundary_vertices
    boundary_normals
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
    depthMap
    depthBoundaries
    refLCOE
    meanWindspeed
    monopileCostModel
end

params = params_struct5(modelset, rotorsamplepointsy, rotorsamplepointsz, turbinez, ambientti,
    rotordiameter, boundaryradius, boundary_vertices, boundary_normals, objectivescale, constraintscaleboundary,
    constraintscalespacing, minimumspacing, hubheight, turbineyaw,
    ctmodels, generatorefficiency, cutinspeed, cutoutspeed, ratedspeed, ratedpower,
    windresource, powermodels, depthMap, depthBoundaries, refLCOE, windspeed, monopileCostModel)


# Set up wrapper functions for the objective and constraints
function boundary_wrapper(x, params)
    # include the relevant params
    boundary_vertices = params.boundary_vertices
    boundary_normals = params.boundary_normals
    constraintscaleboundary = params.constraintscaleboundary

    boundary_x_closed = boundary_vertices[1, :]
    boundary_y_closed = boundary_vertices[2, :]

    boundary_corner_indices = []

    nturbines = Int(length(x)/2)

    # extract the x and y locations of turbines from design variables Vector
    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]

    # get and return boundary distances
    return ff.splined_boundary(turbinex, turbiney, ).*constraintscaleboundary
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
    depthMap = params.depthMap
    refCost = params.refLCOE
    monopileModel = params.monopileCostModel

    # extract x and y locations of turbines from design variables vector
    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]

    xRange = depthBoundaries[1,2] - depthBoundaries[1,1]
    yRange = depthBoundaries[2,2] - depthBoundaries[2,1]

    cost = 0

    for i = eachindex(turbinex)
        i_ind = Int(ceil(10000*turbiney[i]*size(depthMap, 1)))/(yRange*10000)
        up_i_ind = min(size(depthMap, 1), Int(ceil(i_ind)))
        low_i_ind = max(1, Int(floor(i_ind)))

        j_ind = Int(ceil(10000*turbinex[i]*size(depthMap, 2)))/(10000*xRange)
        up_j_ind = min(size(depthMap, 2), Int(ceil(j_ind)))
        low_j_ind = max(1, Int(floor(j_ind)))

        i_indices = [low_i_ind, up_i_ind]
        j_indices = [low_j_ind, up_j_ind]

        depth = depth_at_location(i_ind, j_ind, i_indices, j_indices, depthMap, size(depthMap, 2), size(depthMap, 1))
        
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

    monopilesCost = total_monopile_cost_wrapper(x, params)

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
    #g[1:(end-nturbines)] = spacing_con[:]
    #g[end-nturbines+1:end] = boundary_con[:]
    g[1:(end-(nturbines*6))] = spacing_con[:]
    g[end-(nturbines*6)+1:end] = boundary_con[:]

    obj = coe_wrapper(x, params)[1]

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
monopilesCost = total_monopile_cost_wrapper(currPos, params)
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
best_coe = coe_initial;
best_layout = [copy(turbinex); copy(turbiney)]
best_opt = [copy(turbinex); copy(turbiney)]
# set general lower and upper bounds for design variable
lx = zeros(length(best_layout))
ux = zeros(length(best_layout)) .+ xRange

# set general lower and upper bounds for constraints
ng = Int(nturbines*6 + (nturbines)*(nturbines-1)/2)
lg = [-Inf*ones(Int((nturbines)*(nturbines - 1)/2)); -Inf*ones(nturbines*6)]
ug = [zeros(Int((nturbines)*(nturbines - 1)/2)); zeros(nturbines*6)]

x0 = [copy(turbinex);copy(turbiney)]

# IPOPT options
ip_options = Dict(
    "max_iter" => 50,
    "tol" => 1e-6,
    "print_level" => 0
)
solver = IPOPT(ip_options)

# initialize SNOW options
options = Options(solver=solver, derivatives=ForwardAD())

# -------- RUN OPTIMIZER -----------------------------
# t1 = time()
# xopt
# fopt
# info
# out
# xopt, fopt, info, out = minimize(obj_func!, x0, ng, lx, ux, lg, ug, options)
# t2 = time()
# clk = t2-t1

# coefinal = fopt/objectivescale

# println("Finished in : ", clk, " (s)")
# println("info: ", info)
t1 = time()
for i = 1:10
    global best_coe
    global lx
    global ux
    global ng
    global lg
    global ug
    global best_layout
    global best_opt
    global nturbines
    local turbinex = zeros(nturbines)
    local turbiney = zeros(nturbines)
    local turbinex .+= rand.((0:xRange/2, ))
    local turbiney .+= rand.((0:yRange/2, ))

    local x0 = [copy(turbinex);copy(turbiney)]

    # IPOPT options
    global ip_options
    local solver = IPOPT(ip_options)

    # initialize SNOW options
    local options = Options(solver=solver, derivatives=ForwardAD())

    # -------- RUN OPTIMIZER -----------------------------
    local xopt
    local fopt
    local info
    global out
    t3 = time()
    xopt, fopt, info, out = minimize(obj_func!, x0, ng, lx, ux, lg, ug, options)
    t4 = time()

    local coefinal = fopt/objectivescale
    if coefinal < best_coe
        best_coe = coefinal
        best_layout = x0
        best_opt = xopt
    end
    
    println("TOTAL SOLUTIONS: $i Time: ", t4 - t3)
end
t2 = time()
println()
println()
println("TOTAL TIME: ", t2 - t1)


println("BEST LCOE: $best_coe")
println("BEST LAYOUT: ", best_opt)

# final turbine locations
turbinexopt = copy(best_opt[1:nturbines])
turbineyopt = copy(best_opt[nturbines+1:end])



#-----PLOT OPTIMIZED LAYOUT------------#
fig, ax = plt.subplots(1)

ff.plotlayout!(ax, turbinexopt, turbineyopt, rotordiameter)

ax.set(xlabel="Easting (m)", ylabel="Northing (m)")

#circle = matplotlib.patches.Circle((0.0, 0.0), boundaryradius, fill=false, color="k")
square = matplotlib.patches.Rectangle((0, 0), xRange, yRange, fill=false, color="k")
ax.add_patch(square)

ax.set(xlim=[-xRange*0.01, xRange].*1.01, ylim=[-yRange*0.01, yRange].*1.01)

display(fig)