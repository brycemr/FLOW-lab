using DataFrames
using LinearAlgebra
using Plots
gr()

"""
DownstreamVelocity(u,r0,x)
xr0_ratio - ratio of x(distance downstream) to r0(radius of rotor)
return velocity ratio - returns ratio of wind velocity to ambient air velocity (v/u)

This function applies the Jensen wake model to a single turbine to
calculate the downstream velocity ratio at a given distance. It assumes that
alpha (wake entrainment constant) is equal to 0.1
"""


"""
TEST OUTPUT OF DownstreamVelocity() & DownstreamVelocityRatio()
r = 20 #m
u = 8.1 #m/s
Downwind_Distances = [40, 100] #m
xr0_ratios = [16,10,6]

velocityTable = DataFrame(Downwind_Distances=Int[], Model_Wake_Velocity=Float16[])
velocityRatioTable = DataFrame(x_to_r0_ratio=Float16[], v_to_u_ratio=Float16[])

for distance in Downwind_Distances
    push!(velocityTable, (distance, DownstreamVelocity(u,r,distance)))
end

for ratio in xr0_ratios
    push!(velocityRatioTable, (ratio, DownstreamVelocityRatio(ratio)))
end

println(velocityTable)
println(velocityRatioTable)
"""
function WakeFunction(prevTurbineRatio, xr0_ratio, angle)
    deficit = (1-(prevTurbineRatio/3))*((1/(1+0.1*xr0_ratio))^2)*((1+cosd(9*angle))/2)
    return deficit
end

"""
WindfarmOutputRatio(positions, theta)

positions - an array of all turbine positions in the x-y plane
theta (degrees) - the angle the wind is blowing from the horizontal, i.e. the wind blowing 
in positive x direction would be 0 degrees, positive y would be 90 degrees
rotor_radius - the rotor radius for the turbines used in the wind farm

This function calculates and returns E/E0. this is the total output E for the 
wind farm compared to the output of the same number of non-interacting turbines.

The analysis assumes that the turbines are facing the wind direction
"""
function WindfarmOutputRatio(positions, wind_direction, rotor_radius)
    numTurbines = length(positions)
    #Convert all positions to new coordinate with wind direction as (+) x-axis
    newCoordinates = []
    
    for pos in positions
        newX = pos[1]*cosd(wind_direction) + pos[2]*sind(wind_direction)
        newY = pos[2]*cosd(wind_direction) - pos[1]*sind(wind_direction)
        push!(newCoordinates, [newX, newY])
    end
    
    #Sort the positions in the wind position using a quicksort
    sort!(newCoordinates)

    turbineVelocityRatios = FixedCalculateVelocityRatioAtTurbines(newCoordinates, rotor_radius)[1]
    
    #println("This farm has the following output Ratios: $turbineVelocityRatios")
    output = 0
    for ratio in turbineVelocityRatios
        #output is proportional to wind speed cubed, so ratio^3
        output += ratio^3
    end

    return output/numTurbines
end

"""
Calcuate the Power Output Ratio (compared to a wind farm with 0 interacting turbines) of a wind farm given the positions for the wind farm, rotor radius for the turbines
and wind rose data (an array of speed, frequency, and direction)
"""
function WindfarmYearPerformance(positions, windData, rotor_radius)
    #Power = constants*v^3
    total_output = 0
    possible_output = 0
    for i=axes(windData,1)
        println(windData[i,:])
        windVelocity = windData[i, 1]
        yearlyFraction = windData[i, 2]
        output = (windVelocity^3)*(yearlyFraction)*WindfarmOutputRatio(positions, windData[i, 3], rotor_radius)
        total_output += output
        possible = (windVelocity^3)*(yearlyFraction)
        possible_output += possible
        ratio = output / possible
        println("For this wind data the farm acheives $output of $possible for $ratio")
    end
    finalRatio = total_output/possible_output
    println("In a year we acheive $total_output out of $possible_output for $finalRatio")
    return finalRatio
end

function FindInteractions(upStreamTurbines,currPosition,rotor_radius)
    interactions = []
    #Check backwards so that we later ignore parents starting at nearest turbines
    for i=length(upStreamTurbines):-1:1
        a = (1,0) #horizontal
        b = currPosition - (upStreamTurbines[i] - [(rotor_radius/tand(WAKE_ANGLE)),0]) #This finds the distance from currPosition to the intersection of the top and bottom lines of wake

        angle = acosd(dot(a,b)/(norm(a)*norm(b)))

        if angle <= WAKE_ANGLE
            push!(interactions, i)
        end
    end
    return interactions
end

function VelocityRatioInFarm(x, y, positions, theta, rotor_radius)
    #Convert all positions to new coordinate with wind direction as (+) x-axis
    xs = []
    ys = []
    newCoordinates = []

    calibratedX = x*cosd(theta) + y*sind(theta)
    calibratedY = y*cosd(theta) - x*sind(theta)
    pos = [calibratedX,calibratedY]

    sort!(positions)
    for turbine in positions
        #get x and y's for plotting
        push!(xs, turbine[1])
        push!(ys, turbine[2])
        
        newX = turbine[1]*cosd(theta) - turbine[2]*sind(theta)
        newY = turbine[2]*cosd(theta) + turbine[1]*sind(theta)
        
        if(newX < calibratedX)
            push!(newCoordinates, [newX, newY])
        end
    end

    #Sort the turbines that are in front of the position
    sort!(newCoordinates)
    
    turbineVelocitites, wakeParents = FixedCalculateVelocityRatioAtTurbines(newCoordinates, rotor_radius)
    #println("The velocity ratios in the upstream turbines are $turbineVelocitites")

    velocityRatio = VelocityRatioAtLocation(pos, newCoordinates, turbineVelocitites, wakeParents, rotor_radius)[1]

    PlotWindfarm(x, y, xs, ys, velocityRatio, theta, rotor_radius)

    return velocityRatio
end

function FixedCalculateVelocityRatioAtTurbines(positions,rotor_radius)
    #Velocity Ratio v/u for each turbine in the farm
    Yn = []

    #The first turbine is not in any other turbines wake, so it has a ratio of 1
    push!(Yn, 1)
    numTurbines = length(positions)
    wakeParents = [[0]]

    for i=2:numTurbines
        velocityRatio, newParents = VelocityRatioAtLocation(positions[i], positions[1:i-1], Yn, wakeParents, rotor_radius)

        push!(wakeParents, newParents)
        push!(Yn, velocityRatio)     
    end
    return Yn, wakeParents
end

function VelocityRatioAtLocation(location, upstreamTurbines, upstreamVelocities, wakeParents, rotor_radius)
    interactions = FindInteractions(upstreamTurbines, location, rotor_radius)

    velocityRatio = 1
    
    ignoreParents=Set()
    ##COMMENTS
    for j in interactions
        for parent in wakeParents[j]
            push!(ignoreParents, parent)
        end

        if !(j in ignoreParents)
            xr0_ratio = norm(location - upstreamTurbines[j]) / rotor_radius
            velocityRatio -= WakeFunction(upstreamVelocities[j], xr0_ratio, 0) 
        end
    end
    return velocityRatio, interactions
end

function PlotWindfarm(x, y, xs, ys, velocityRatio, theta, rotor_radius)
    scatter(xs, ys, seriestype=:scatter, title = "Wind Farm Turbine Layout", label="", markershape=:cross, size=(400,400), markerstrokewidth=4)
    outputRatio = round(velocityRatio, digits=3)

    for i=eachindex(xs)
        turbineTop = [xs[i]+rotor_radius*cosd(theta+90), ys[i]+rotor_radius*sind(theta+90)]
        turbineBottom = [xs[i]+rotor_radius*cosd(theta-90), ys[i]+rotor_radius*sind(theta-90)]
        wakeTop = [turbineTop[1]+100*rotor_radius*cosd(theta), turbineTop[2]+100*rotor_radius*tand(theta+WAKE_ANGLE)] #100*rotor_radius because at x = 100r0 the velocity is more than .99u
        wakeBottom = [turbineBottom[1]+100*rotor_radius*cosd(theta), turbineBottom[2]+100*rotor_radius*tand(theta-WAKE_ANGLE)] 
        wakeXs = [wakeBottom[1], turbineBottom[1], turbineTop[1], wakeTop[1]]
        wakeYs = [wakeBottom[2], turbineBottom[2], turbineTop[2], wakeTop[2]]
        plot!(wakeXs, wakeYs, linealpha=0.50, label="")
    end

    scatter!([x],[y], color="red", markershape=:xcross, label="$outputRatio u")
    tail = [x-50*cosd(theta), x]
    head = [y-50*sind(theta), y]
    plot!(tail, head, linealpha=0.50, color="blue", label="", arrow=true)
    plot!(xlims=(-110,110), framestyle=:none)
    ylims!((-110,110))
    
    filename = "Windfarm_Plot.pdf"
    savefig(filename)
end


#INPUTS TO SOLVE FOR
WAKE_ANGLE = 20

#Create Circle
circle = []
circle_radius = 100
numTurbines = 10
angle = 360/numTurbines
for i=1:numTurbines
    push!(circle,[circle_radius*cosd((angle*i)+10), circle_radius*sind((angle*i)+10)])
end

#Create straight line of turbines
line=[[0,0],[50,0],[100,0],[150,0],[200,0],[250,0],[300,0],[350,0],[400,0],[450,0]]


#Create Wind Rose (velocity, fequency, direction)
windRose = [
    14 0.1 0;
    10 0.3 90;
    5 0.4 180;
    20 0.2 270]

#Call Functions
#println("Wind farm output: ")
#WindfarmYearPerformance(line, windRose, 10)


#TEST RANGE OF ANGLES FOR CIRCLE
#for i=-18:18
    #print("For angle $i the output ratio is: ")
    #println(WindfarmOutputRatio(circle, i, 10))
#end


output = VelocityRatioInFarm(0, 0, circle, 30, 10)
#println("The output of this farm is $output")
