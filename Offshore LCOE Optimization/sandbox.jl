using Plots
using Distributions
include("Monopile.jl")
include("Depths.jl")

println("Start")
import_depth_from_xyz("C:\\Users\\bryce\\Desktop\\FLOW-lab\\Offshore LCOE Optimization\\H06155.xyz", 100, 100, [-90.4288702*1.0001, -90.3360071*.9999], [28.6673358*0.9999, 28.7491368*1.0001], "linear")
println("Finished")
