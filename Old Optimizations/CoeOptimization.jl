using FLOWFarm; const ff = FLOWFarm
using PyPlot; const plt = PyPlot
using VectorizedRoutines.Matlab: meshgrid
using SNOW
using Plots

turbinex = [-500.0, -500.0, -500.0, 0.0, 0.0, 0.0, 500.0, 500.0, 500.0]
turbiney = [-500.0, 0.0, 500.0, -500.0, 0.0, 500.0, -500.0, 0.0, 500.0]
#turbiney = [0.0, 0.0, 0.0]
#turbinex = [-240.0, 0, 240]
#turbinex = [-500.0, -500.0, -500.0, -500.0, -500.0, -250.0, -250.0, -250.0, -250.0, -250.0, 0.0, 0.0, 0.0, 0.0, 0.0, 250.0, 250.0, 250.0, 250.0, 250.0, 500.0, 500.0, 500.0, 500.0, 500.0]
#turbiney = [-500.0, -250.0, 0.0, 250.0, 500.0, -500.0, -250.0, 0.0, 250.0, 500.0, -500.0, -250.0, 0.0, 250.0, 500.0, -500.0, -250.0, 0.0, 250.0, 500.0, -500.0, -250.0, 0.0, 250.0, 500.0]


# INITIAL TURBINE LOCATIONS

nturbines = length(turbinex)

turbinez = zeros(nturbines)

turbineyaw = zeros(nturbines)

# set wind farm boundary parameters in meters
boundarycenter = [0.0, 0.0]
boundaryradius = 900

boundary_vertices = [-boundaryradius -boundaryradius; boundaryradius -boundaryradius; boundaryradius boundaryradius; -boundaryradius boundaryradius]
boundary_normals = [-1 0; 0 -1; 1 0; 0 1]

# set wind turbine design parameters
rotordiameter = zeros(nturbines) .+ 125 # m 
hubheight = zeros(nturbines) .+ 90      # m
cutinspeed = zeros(nturbines) .+ 3.0    # m/s
cutoutspeed = zeros(nturbines) .+ 25.0   # m/s 
ratedspeed = zeros(nturbines) .+ 11.4   # m/s 
ratedpower = zeros(nturbines) .+ 5.0E6  # W (5 MW)
generatorefficiency = ones(nturbines)

# VISUALIZING THE WIND FARM LAYOUT


fig, ax = plt.subplots(1)

ff.plotlayout!(ax, turbinex, turbiney, rotordiameter)

ax.set(xlabel="Easting (m)", ylabel="Northing (m)")

#circle = matplotlib.patches.Circle((0.0, 0.0), boundaryradius, fill=false, color="k")
square = matplotlib.patches.Rectangle((-1000.0, -1000.0), 2000, 2000, fill=false, color="k")
ax.add_patch(square)

ax.set(xlim=[-boundaryradius, boundaryradius].*1.01, ylim=[-boundaryradius, boundaryradius].*1.01)

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
display(ff.plotwindresource!(windresource))


# SETTING UP ANALYSIS MODELS
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

# RUNNING THE ANALYSIS

# Annual Energy production aep
aep = ff.calculate_aep(turbinex, turbiney, turbinez, rotordiameter,
    hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
    cutoutspeed, ratedspeed, ratedpower, windresource, powermodels,
    modelset, rotor_sample_points_y=rotorsamplepointsy, 
    rotor_sample_points_z=rotorsamplepointsz)


refTCC = 1952.0
refBOS = 4420 + 474.0
refFC = 672.0  # $/kW
refFCR = 0.0764 # /year
refCapEx = 140.56
cost = ff.Levelized(refTCC, refBOS, refFC, refFCR, refCapEx)
total_ratedpower = sum(ratedpower)
coe_aep = aep/total_ratedpower # For the LCOE equation we need aep in MWh/MW/year not Wh/year (same as hr/year)
coe_ratedpower = ratedpower./1000 # needs to be in units of kw
coe = ff.cost_of_energy(rotordiameter, hubheight, coe_ratedpower, coe_aep, cost)

println("$aep Watt-hours per year")
# AEP in each direction
state_aeps = ff.calculate_state_aeps(turbinex, turbiney, turbinez, rotordiameter,
        hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
        cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
        rotor_sample_points_y=rotorsamplepointsy, rotor_sample_points_z=rotorsamplepointsz,
        hours_per_year=365.25*24.0, weighted=true)

println("$state_aeps Watt-hours per year")

# --------SETTING UP CONSTRAINTS AND OPTIMIZATION FUNCTIONS-------------

# scale objective derivatives to be between 0 and 1
objectivescale = 1E5

# scale boundary constraint derivatives to be between 0 and 1
constraintscaleboundary = 1.0E-3

# scale spacing constraint derivatives to be between 0 and 1
constraintscalespacing = 1.0

# set the minimum spacing between turbines 
minimumspacing = 125*1.2

# set up a struct for use in optimization functions, these are the non-differentiated parameters
mutable struct params_struct4{}
    modelset
    rotorsamplepointsy
    rotorsamplepointsz
    turbinez
    ambientti
    rotordiameter
    boundarycenter
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
end

params = params_struct4(modelset, rotorsamplepointsy, rotorsamplepointsz, turbinez, ambientti,
    rotordiameter, boundarycenter, boundaryradius, boundary_vertices, boundary_normals, objectivescale, constraintscaleboundary,
    constraintscalespacing, minimumspacing, hubheight, turbineyaw,
    ctmodels, generatorefficiency, cutinspeed, cutoutspeed, ratedspeed, ratedpower,
    windresource, powermodels)

# Set up wrapper functions for the objective and constraints
function boundary_wrapper(x, params)
    # include the relevant params
    boundarycenter = params.boundarycenter
    boundaryradius = params.boundaryradius
    boundary_vertices = params.boundary_vertices
    boundary_normals = params.boundary_normals
    constraintscaleboundary = params.constraintscaleboundary

    nturbines = Int(length(x)/2)

    # extract the x and y locations of turbines from design variables Vector
    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]

    # get and return boundary distances
    #return ff.circle_boundary(boundarycenter, boundaryradius, turbinex, turbiney).*constraintscaleboundary
    return ff.convex_boundary(boundary_vertices, boundary_normals, turbinex, turbiney).*constraintscaleboundary
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

    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]

    # calculate_aep
    aep = ff.calculate_aep(turbinex, turbiney, turbinez, rotordiameter,
    hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
    cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
    rotor_sample_points_y=rotorsamplepointsy,rotor_sample_points_z=rotorsamplepointsz)

    total_ratedpower = sum(ratedpower)
    coe_aep = aep/total_ratedpower # For the LCOE equation we need aep in MWh/MW/year not Wh/year (same as hr/year)
    coe_ratedpower = ratedpower./1000 # needs to be in units of kw
    
    refTCC = 1952.0
    refBOS = 4420 + 474.0
    refFC = 672.0  # $/kW
    refFCR = 0.0764 # /year
    refCapEx = 140.56
    
    cost = ff.Levelized(refTCC, refBOS, refFC, refFCR, refCapEx)
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
    g[1:(end-(nturbines*4))] = spacing_con[:]
    g[end-(nturbines*4)+1:end] = boundary_con[:]

    obj = coe_wrapper(x, params)[1]
    
    return obj
end

# generate objective function wrapper
obj_func!(g, x) = wind_farm_opt!(g, x, params)


# ------------SETUP THE OPTIMIZER --------------------------------------

# initialize design variable vector
x0 = [copy(turbinex);copy(turbiney)]

# set general lower and upper bounds for design variable
lx = zeros(length(x0)) .- boundaryradius
ux = zeros(length(x0)) .+ boundaryradius

# set general lower and upper bounds for constraints
ng = Int(nturbines*4 + (nturbines)*(nturbines-1)/2)
lg = [-Inf*ones(Int((nturbines)*(nturbines - 1)/2)); -Inf*ones(nturbines*4)]
ug = [zeros(Int((nturbines)*(nturbines - 1)/2)); zeros(nturbines*4)]

# IPOPT options
ip_options = Dict(
    "max_iter" => 100,
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
println("Initial COE: ", coe)
println("Final COE: ", coefinal)
println("COE improvement (%) = ", -100*(coefinal - coe)/coe)
println("Final Turbine Positioning (x, y):")

# final turbine locations
turbinexopt = copy(xopt[1:nturbines])
turbineyopt = copy(xopt[nturbines+1:end])
println(turbinexopt)
println(turbineyopt)


#-----OPTIMIZED LAYOUT------------
fig, ax = plt.subplots(1)

ff.plotlayout!(ax, turbinexopt, turbineyopt, rotordiameter)

ax.set(xlabel="Easting (m)", ylabel="Northing (m)")

# and the wind farm boundary
square = matplotlib.patches.Rectangle((-1000.0, -1000.0), 2000, 2000, fill=false, color="k")
ax.add_patch(square)

# set limits on the plot region
ax.set(xlim=[-boundaryradius, boundaryradius].*1.01, ylim=[-boundaryradius, boundaryradius].*1.01)

display(fig)
