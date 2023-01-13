"""
    mapDepths(resolution, x_boundaries, y_boundaries, distribution, gausStartDepth, linnStartDepth, linDepthVariance, centerIsShallow, isLinear, isGaus1D, isGaus2D, gausMaxDepth, gausMinDepth)

Creates a heat map for a specified water depth pattern.
# Arguments
-'resolution::Int': Resolution of the of the x and y axis for the heat map
-'x_boundaries::array': Vector of the minimum and maximum of the x-axis
-'y_boundaries::array': Vector of the minimum and maximum of the y-axis
-'distribution::??': ?????
-
"""
function mapDepths(resolution, x_boundaries, y_boundaries, distribution, gausStartDepth, linStartDepth, linDepthVariance, centerIsShallow, isLinear, isGaus1D, isGaus2D, gausMaxDepth, gausMinDepth)
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
            index = i*dimY + j + 1
            pos = [lowerX+i*xstep, lowerY+j*ystep]

            if isGaus2D
                curr_depth = gaussianTwoD(distribution, gausMaxDepth, gausMinDepth, pos, centerIsShallow)
            else
                linearQuant = (i+0.001)/dimX
                gaussQuant = (j+0.001)/dimY

                distDepth = quantile(distribution, gaussQuant)[1]
                difference = abs(gausStartDepth - distDepth)
                curr_depth = (isGaus1D ? gaussianOneD(gausStartDepth, difference, centerIsShallow) : 0)
                curr_depth = (isLinear ? curr_depth + linearDepth(linDepthVariance, linStartDepth, linearQuant) : curr_depth)
            end

            depths[:,index] = [pos[1]; pos[2]; curr_depth]
        end
    end

    xs = depths[1,:]
    ys = depths[2,:]
    zs = depths[3,:]

    x, y = sort(unique(xs)), sort(unique(ys)) 

    C = fill(NaN, length(x), length(y)) 
    for i = eachindex(x)
        for j = eachindex(y)
            C[i,j] = zs[(i-1)*length(y)+j]
        end
    end

    heat_map = heatmap(x, y, C, clims=(10,60)) 
    display(heat_map)
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
    if centerIsShallow
        depth = maxDepth - depth_diff*percentile/meanPercentile 
    else
        depth = minDepth + depth_diff*percentile/meanPercentile
    end

    return depth
end