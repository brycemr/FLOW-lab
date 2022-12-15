using FLOWFarm; const ff = FLOWFarm
using PyPlot; const plt = PyPlot
using VectorizedRoutines.Matlab: meshgrid
using SNOW

turbinex = [-240.0, -240.0, -240.0, 0.0, 0.0, 0.0, 240.0, 240.0, 240.0]
turbiney = [-240.0, 0.0, 240.0, -240.0, 0.0, 240.0, -240.0, 0.0, 240.0]

nturbines = length(turbinex)

turbinez = zeros(nturbines)

turbineyaw = zeros(nturbines)

# set wind farm boundary parameters in meters
boundarycenter = [0.0, 0.0]
boundaryradius = hypot(300,300)

# set wind turbine design parameters
rotordiameter = zeros(nturbines) .+ 80.0 # m 
hubheight = zeros(nturbines) .+ 70      # m
cutinspeed = zeros(nturbines) .+ 4.0    # m/s
cutoutspeed = zeros(nturbines) .+ 25.0   # m/s 
ratedspeed = zeros(nturbines) .+ 16.0   # m/s 
ratedpower = zeros(nturbines) .+ 2.0E6  # W
generatorefficiency = ones(nturbines)

# VISUALIZING THE WIND FARM LAYOUT

fig, ax = plt.subplots(1)

ff.plotlayout!(ax, turbinex, turbiney, rotordiameter)

ax.set(xlabel="Easting (m)", ylabel="Northing (m)")

#circle = matplotlib.patches.Circle((0.0, 0.0), boundaryradius, fill=false, color="k")
#ax.add_patch(circle)

ax.set(xlim=[-boundaryradius, boundaryradius].*1.01, ylim=[-boundaryradius, boundaryradius].*1.01)

display(fig)

# get sample points for rotor swept area for determining inflow wind speed
nsamplepoints = 50
rotorsamplepointsy, rotorsamplepointsz = ff.rotor_sample_points(nsamplepoints, method="sunflower")

# SETUP WIND RESOURCE
windspeed = 8.0 # m/s 
airdensity = 1.1716 # kg/m^3
ambientti = 0.1 # %
shearexponent = 0.15
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
objectivescale = 1E-6

# scale boundary constraint derivatives to be between 0 and 1
constraintscaleboundary = 1.0E-3

# scale spacing constraint derivatives to be between 0 and 1
constraintscalespacing = 1.0

# set the minimum spacing between turbines 
minimumspacing = 100

# set up a struct for use in optimization functions, these are the non-differentiated parameters
mutable struct params_struct2{}
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
end

params = params_struct2(modelset, rotorsamplepointsy, rotorsamplepointsz, turbinez, ambientti,
    rotordiameter, boundarycenter, boundaryradius, objectivescale, constraintscaleboundary,
    constraintscalespacing, minimumspacing, hubheight, turbineyaw,
    ctmodels, generatorefficiency, cutinspeed, cutoutspeed, ratedspeed, ratedpower,
    windresource, powermodels)

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

# Set up the aep wrapper function
function aep_wrapper(x, params)
    
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
    aep = objectivescale*ff.calculate_aep(turbinex, turbiney, turbinez, rotordiameter,
    hubheight, turbineyaw, ctmodels, generatorefficiency, cutinspeed,
    cutoutspeed, ratedspeed, ratedpower, windresource, powermodels, modelset,
    rotor_sample_points_y=rotorsamplepointsy,rotor_sample_points_z=rotorsamplepointsz)

    return aep
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

    obj = -aep_wrapper(x, params)[1]
    
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

aepfinal = -fopt/objectivescale

println("Finished in : ", clk, " (s)")
println("info: ", info)
println("Initial AEP: ", aep)
println("Final AEP: ", aepfinal)
println("AEP improvement (%) = ", 100*(aepfinal - aep)/aep)

# final turbine locations
turbinexopt = copy(xopt[1:nturbines])
turbineyopt = copy(xopt[nturbines+1:end])


#-----OPTIMIZED LAYOUT------------
fig, ax = plt.subplots(1)

ff.plotlayout!(ax, turbinexopt, turbineyopt, rotordiameter)

ax.set(xlabel="Easting (m)", ylabel="Northing(m)")

# and the wind farm boundary
circle = matplotlib.patches.Circle((0.0, 0.0), boundaryradius, fill=false, color="k")
ax.add_patch(circle)

# set limits on the plot region
ax.set(xlim=[-boundaryradius, boundaryradius].*1.01, ylim=[-boundaryradius, boundaryradius].*1.01)

display(fig)
