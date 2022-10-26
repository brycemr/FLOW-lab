# Store Bathymetry data as x,y,z array, normalize it so that each index corresponds to 
# a set step in x and y (average distance) O(xy) or must be provided in this format, in which 
# case O(1) to caclulate step size, max and min
data = zeros(10,10)
depth = 15
for i = 1:10
    for j = 1:10
        data[i,j] = depth 
    end
    global depth = depth + 10
end

println(data)
# while turbine is being positioned round x and y to nearest square (divide by step) (O(1))
turbinex = 93
turbiney = 65
step = 10
x_index = floor(Int, turbinex / step)
y_index = floor(Int, turbiney / step)
println(x_index)
println(y_index)
# create a lookup table, input depth values and cost values, (min_depth, max_depth, cost) for each turbine
# O(t) t is num turbines
struct Foundation
    name
    min_depth
    max_depth
    cost
end
onshore = Foundation("onshore", 10, 60, 100)
offshore = Foundation("offshore", 40, 200, 150)
turbines = [onshore, offshore]
# for z value at given square, lookup associated cost for required turbine, 
# if multiple turbines are possible, choose lowest cost, add to calculation of LCOE or COE
minCost = 200 
foundationType = "none"
for t in turbines
    if t.max_depth > data[x_index,y_index] > t.min_depth
        if t.cost < minCost
            global minCost = cost
            global foundationType = t.name
        end
    end
end

println(foundationType)
println(minCost)

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

function monopile_cost_wrapper(x, params)
    cost = 0
    
    nturbines = Int(length(x)/2)
    depthMap = params.depthMap

    # extract x and y locations of turbines from design variables vector
    turbinex = x[1:nturbines]
    turbiney = x[nturbines+1:end]
    
    for i = eachindex(turbinex)
        # This works for finding the x and y index of a rectangular depth map with constant x and y steps
        depth_x_index = ceil((turbinex[i] - minDepthMapX) / xStepSize) 
        depth_y_index = ceil((turbiney[i] - minDepthMapY) / yStepSize)
        depth = depthMap[depth_x_index, depth_y_index]
        cost = cost + calculate_monopile_cost(depth, params)
    end

    return cost
end

function calculate_monopile_cost(depth, params)
    # place holder monopile cost equation 600 tons / 25 meters => mass ~= (24*depth) tons  
    # steel is ~$700 per ton, this is just material does not consider 
    # manufacturing costs, transportation costs, etc.
    # based on figure 3 of NREL monopile vs jacket paper
    return (depth*24)*700
end


function coe_wrapper(x, params)
    cost = monopile_cost_wrapper(x, params)
    aep = aep_wrapper(x, params)
    coe = cost/aep
    return coe
end

