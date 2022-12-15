#=
SAMPLE HEATMAP

using Plots
a = [1, 2, 3, 1, 2, 3, 1, 2]
b = [1, 1, 1, 2, 2, 2, 3, 3]
c = [1, 5, 4, 3, 4, 2, 1, 3]
x, y = sort(unique(a)), sort(unique(b)) # the lattice x and y
C = fill(NaN, length(y), length(x)) # make a 2D C with NaNs
for (i,j,v) in zip(a,b,c) # fill C with c values
    C[j,i] = v
end
println(x)
println(y)
println(C)
heatmap(x, y, C, clims=(0,5)) # success!
=#
#import Pkg
#Pkg.add("Roots")
using Distributions
using Roots

println(quantile(data, [0.98]))

function mapDepths(resolution, x_boundaries, y_boundaries, distribution, gausStartDepth, linStartDepth, linDepthVariance, centerIsShallow, gaus2D, gausMaxDepth, gausMinDepth)
    upperX = x_boundaries[2]
    lowerX = x_boundaries[1]
    upperY = y_boundaries[2]
    lowerY = y_boundaries[1]
    
    xstep = (upperX - lowerX) / resolution
    ystep = (upperY - lowerY) / resolution
    dimX = Int((upperX - lowerX) / xstep)+1
    dimY = Int((upperY - lowerY) / ystep)+1
    depths = Array{Float64}(undef, 3, (dimX)*(dimY))
    println("xstep: $xstep, ystep: $ystep, dimX: $dimX, dimY: $dimY")
    
    for i = 0:dimX-1
        for j = 0:dimY-1
            # Currently generates random depth data, instead would find depth data at this point
            index = i*dimY + j + 1
            # COULD ADD INPUTS THAT ALLOW SPECIFICATION OF LINEAR ON X OR Y AND SAME FOR GAUSS
            pos = [lowerX+i*xstep, lowerY+j*ystep]

            if gaus2D
                curr_depth = gaussianTwoD(distribution, gausMaxDepth, gausMinDepth, pos, centerIsShallow)
                print("pos: $pos, curr_depth: $curr_depth ")
            else
                linearQuant = (i+0.001)/dimX
                gaussQuant = (j+0.001)/dimY

                distDepth = quantile(distribution, gaussQuant)[1]
                difference = abs(gausStartDepth - distDepth)
                curr_depth = (gausDepth ? gaussianOneD(gausStartDepth, difference, centerIsShallow) : 0)
                curr_depth = (linDepth ? curr_depth + linearDepth(linDepthVariance, linStartDepth, linearQuant) : curr_depth)
            end

            depths[:,index] = [pos[1]; pos[2]; curr_depth]
        end
    end
    
    return depths
end

function linearDepth(depthVariance, baseDepth, quantile)
    return baseDepth + depthVariance*quantile
end

function gaussianOneD(data_mean, difference, centerIsShallow)
    if centerIsShallow
        return data_mean + difference
    else
        return data_mean - difference
    end
end

function gaussianTwoD(distribution, maxDepth, minDepth, pos, centerIsShallow)
    depth_diff = maxDepth - minDepth
    meanPercentile = pdf(distribution, [0.0, 0.0])
    percentile = pdf(distribution, pos)
    quantile = percentile/meanPercentile
    println(" depth_diff: $depth_diff, quantile: $quantile")
    if centerIsShallow
        depth = maxDepth - depth_diff*percentile/meanPercentile 
    else
        depth = minDepth + depth_diff*percentile/meanPercentile
    end

    return depth
end

gaus1D = true
gaus2D = true
gausDepthOutward = true
linDepth = false
gausMaxDepth = 60
gausMinDepth = 10
linStartDepth = 0
linDepthVariance = 0

depthBoundaries = [-1250.0 1250.0; -1250.0 1250.0]

if gaus2D
    means = [0.0, 0.0] # Use center x, and y pos as means
    coVariants = [5.0e5 4.0e5; 4.0e5 5.0e5] # 
    distribution = MvNormal(means, coVariants)
else
    gausStartDepth = 10
    depthStdDev = (gausMaxDepth - gausMinDepth)/3
    distribution = Normal(gausStartDepth, depthStdDev)
end



#depths = mapDepths(60, [-1050, 1050], [-1050, 1050], distribution, gausStartDepth, linStartDepth, linDepthVariance, gausDepthOutward, gaus2D, gausMaxDepth, gausMinDepth)
# Plot Depth Data
xs = depths[1,:]
ys = depths[2,:]
zs = depths[3,:]

x, y = sort(unique(xs)), sort(unique(ys)) # the lattice x and y

println()
#println(x)
println()
#println(y)
#println(zs)

C = fill(NaN, length(x), length(y)) # make a 2D C with NaNs
for i = eachindex(x)
    for j = eachindex(y)
        # Currently generates random depth data, instead would find depth data at this point
        C[i,j] = zs[(i-1)*length(y)+j]
    end
end
println()
#println(C)
#heatmap(x, y, C, clims=(10,30)) # success!

#x = [i/100 for i in 0:100]
#y = [j/100 for j in 0:100]

#Z = [pdf(d,[i/100,j/100]) for i in 0:100, j in 0:100]
#print(Z)
#heatmap(x, y, C, clims=(0,60))

#m = [1.0, 1.0]
#C = [0.1 0; 0 0.1]
#d = MvNormal(m, C)
#println(gaussianTwoD(d, 40, 10, [0.5, 0.5]))
#println(pdf(d, [0.5, 0.0]))
#println(pdf(d, [0.0, 0.5]))
#println(pdf(d, [0.5, 0.5]))
#println(pdf(d, [1.0, 0]))
#println(pdf(d, [1.0, 1.0]))
#gp2 = GP(X,y,m,se)

#p1 = plot(gp2; title="Mean of GP")
#p2 = plot(gp2; var=true, title="Variance of GP", fill = true)
#plot(p1, p2, fmt=:png)

# CODE ORGANIZATION
# Create the gaussian Distributions once, create a function that calculates the depth based on max depth, min depth, distribution, and curr pos
function monopile_diam_equation(Dp, data)
    # Equation 99 and 101
    yield_stress = data.yield_stress
    material_factor = data.material_factor 
    M_50y = data.M_50y
    A = (yield_stress * pi) / (4 * material_factor * M_50y)
    res = A * ((0.99 * Dp - 0.00635) ^ 3) * (0.00635 + 0.01 * Dp) - Dp

    return res
end
yield_stress = 30000000
material_factor =  1.1
M_50y = 1000000
A = (yield_stress * pi) / (4 * material_factor * M_50y)
    
g(Dp) = res = A * ((0.99 * Dp - 0.00635) ^ 3) * (0.00635 + 0.01 * Dp) - Dp
f(x) = x^5 - x + 1/2
sol = find_zero(g, 10)
println(sol)