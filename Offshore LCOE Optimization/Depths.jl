
function mapDepths(resolution, x_boundaries, y_boundaries, distribution, gausStartDepth, linStartDepth, linDepthVariance, centerIsShallow, isLinear, isGaus1D, isGaus2D, gausMaxDepth, gausMinDepth)
    upperX = x_boundaries[2]
    lowerX = x_boundaries[1]
    upperY = y_boundaries[2]
    lowerY = y_boundaries[1]
    
    xstep = (upperX - lowerX) / resolution
    ystep = (upperY - lowerY) / resolution
    dimX = Int(round((upperX - lowerX) / xstep))+1
    dimY = Int(round((upperY - lowerY) / ystep))+1
    depths = Array{Float64}(undef, 3, (dimX)*(dimY))
    
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


# At beginning of program latitude and longitude should be mapped relative to center of the wind farm (0, 0) to form the x and y coordinates. A grid of depths can then 
# be created with these x and y coordinates, if a coordinate does not have a depth provided it can be averages from neighbors in x and y (average of all four neighbors)
# when looking up depths, interpolate between points for x and y
# Take in depth data that has been split into three different arrays
function map_depth_data(depth_x_range, depth_y_range, depths, length_x, length_y)
    depths = Array{Float64}(undef, 3, (dimX)*(dimY))
end


# May be able to find indices by doing integer division with resolution on distance from min edge
function depth_at_location(i_indx, j_indx, i_indices, j_indices, depth_map, resolution_x, resolution_y)
    corner_1_depth = depth_map[i_indices[1], j_indices[1]]
    corner_2_depth = depth_map[i_indices[2], j_indices[1]]
    corner_3_depth = depth_map[i_indices[1], j_indices[2]]
    corner_4_depth = depth_map[i_indices[2], j_indices[2]]

    first_x_interpolation = interpolate_depth(j_indx - j_indices[1], corner_1_depth, corner_2_depth)
    second_x_interpolation = interpolate_depth(j_indx - j_indices[1], corner_3_depth, corner_4_depth)
    
    depth = interpolate_depth(i_indx-i_indices[1], first_x_interpolation, second_x_interpolation)
    
    return depth
end

function interpolate_depth(x, depth1, depth2)
    return depth1 + ((depth2 - depth1)) * x
end

# Create a 2D array that stores each depth entry with the latitude, longitude, and depth. Plot a heat_map
# Resolution specifies how many cells in x and y
function import_depth_from_xyz(filepath, num_cells_x, num_cells_y, min_max_latitude, min_max_longitude, fill_kernel="gaussian")
    f = open(filepath, "r")

    line = readline(f)
    header = split(line, ",")
    println(line)
    println(header)
    println(size(header))

    # Locate the columns that hold lat, long, and depth
    lat_col = 0
    long_col = 0
    depth_col = 0
    min_depth = Inf
    max_depth = 0

    for i in eachindex(header)
        if header[i] == "lat"
            lat_col = i
        elseif header[i] == "long"
            long_col = i
        elseif header[i] == "depth"
            depth_col = i
        end
    end
    println(lat_col)
    println(long_col)
    println(depth_col)

    # Identify the maximum and minimum lat, long, and depth. 
    min_lat = min_max_latitude[1]
    max_lat = min_max_latitude[2]
    min_long = min_max_longitude[1]
    max_long = min_max_longitude[2]
    
    # Get the range in meters
    N_S_range = latitude_to_meters(max_lat-min_lat)
    E_W_range = longitude_to_meters((max_lat+min_lat)/2, max_long-min_long)

    println(max_lat," ", min_lat)
    println(max_long," ", min_long)
    println(N_S_range)
    println(E_W_range)

    # Create x and y values for heat map
    xs = zeros(num_cells_x)
    ys = zeros(num_cells_y)
    for i = 1:num_cells_x
        xs[i] = i*E_W_range/num_cells_x
    end
    for j = 1:num_cells_y
        ys[j] = j*N_S_range/num_cells_y
    end
    
    # Create 2D array with specified num_cells for latitude and longitude
    depth_map = fill(0.0, num_cells_y, num_cells_x) 
    num_entries = fill(0, num_cells_y, num_cells_x)
    # Parse all the data points again, filling in depth array
    for line in readlines(f)
        data = split(line, ",")
        lat = parse(Float64, data[lat_col])
        long = parse(Float64, data[long_col])
        depth = parse(Float64, data[depth_col])

        if lat > min_lat && lat < max_lat && long > min_long && long < max_long
            if depth > max_depth
                max_depth = depth
            end
            if depth < min_depth
                min_depth = depth
            end
            
            # Convert to m offset from bottom left corner (SW corner), and turn this into index
            x = longitude_to_meters(lat, long - min_long)
            y = latitude_to_meters(lat - min_lat)

            x_ind = min(num_cells_x, Int(round(x * num_cells_x / E_W_range))+1)
            y_ind = min(num_cells_y, Int(round(y * num_cells_y / N_S_range))+1)

            # If more than one depth entry falls into a cell, take average of the values
            n_entry = num_entries[y_ind, x_ind]
            num_entries[y_ind, x_ind] = n_entry + 1
            prev_depth = depth_map[y_ind, x_ind]
            depth_map[y_ind, x_ind] = (depth + prev_depth*n_entry)/(n_entry + 1)
        end
    end

    # Display initial depth map
    heat_map = heatmap(xs, ys, depth_map, clims=(round(min_depth*0.75), round(max_depth*1.25))) 
    display(heat_map)

    # Fill in zero cells
    if fill_kernel == "linear"
        filled_map = linear_map_fill(depth_map, num_cells_x, num_cells_y)
    else
        filled_map = gaussian_map_fill(depth_map, num_cells_x, num_cells_y)
    end

    heat_map = heatmap(xs, ys, filled_map, clims=(round(min_depth*0.75), round(max_depth*1.25))) 
    display(heat_map)

    return filled_map
end

function latitude_to_meters(diff_latitude)
    return diff_latitude*111.32*1000
end

# Take in latitude as degrees
function longitude_to_meters(latitude, diff_longitude)
    return diff_longitude*40075*1000*cos(latitude * pi / 180) / 360
end

function gaussian_map_fill(depth_map, num_cells_x, num_cells_y)
    new_map = deepcopy(depth_map)
    repeat = true
    while repeat
        println("ITERATION")
        repeat = false
        for i = 1:num_cells_y
            for j = 1:num_cells_x
                sum = 0.0
                count = 0
                if new_map[i, j] == 0.0
                    if i > 1
                        sum += 2 * new_map[i-1, j]
                        if new_map[i-1, j] > 0
                            count += 2
                        end

                        # Blur with upper left and lower left cells
                        if j > 1
                            sum += new_map[i, j-1]
                            if new_map[i, j-1] > 0
                                count += 1
                            end
                        end
                        if j < num_cells_x
                            sum += new_map[i, j+1]
                            if new_map[i, j+1] > 0
                                count += 1
                            end
                        end

                    end
                    if i < num_cells_y
                        sum += 2 * new_map[i+1, j]
                        if new_map[i+1, j] > 0
                            count += 2
                        end

                        # Blur with upper right and lower right cells
                        if j > 1
                            sum += new_map[i, j-1]
                            if new_map[i, j-1] > 0
                                count += 1
                            end
                        end
                        if j < num_cells_x
                            sum += new_map[i, j+1]
                            if new_map[i, j+1] > 0
                                count += 1
                            end
                        end
                    end
                    if j > 1
                        sum += 2 * new_map[i, j-1]
                        if new_map[i, j-1] > 0
                            count += 2
                        end
                    end
                    if j < num_cells_x
                        sum += 2 * new_map[i, j+1]
                        if new_map[i, j+1] > 0
                            count += 2
                        end
                    end

                    if count == 0
                        repeat = true
                    else
                        new_map[i, j] = sum/count
                    end
                    
                end
            end
        end
    end
    return new_map
end

function linear_map_fill(depth_map, num_cells_x, num_cells_y)
    new_map = deepcopy(depth_map)
    repeat = true
    while repeat
        println("ITERATION")
        repeat = false
        for i = 1:num_cells_y
            for j = 1:num_cells_x
                sum = 0.0
                count = 0
                if new_map[i, j] == 0.0
                    if i > 1
                        sum += new_map[i-1, j]
                        if new_map[i-1, j] > 0
                            count += 1
                        end
                    end
                    if i < num_cells_y
                        sum += new_map[i+1, j]
                        if new_map[i+1, j] > 0
                            count += 1
                        end
                    end
                    if j > 1
                        sum += new_map[i, j-1]
                        if new_map[i, j-1] > 0
                            count += 1
                        end
                    end
                    if j < num_cells_x
                        sum += new_map[i, j+1]
                        if new_map[i, j+1] > 0
                            count += 1
                        end
                    end

                    if count == 0
                        repeat = true
                    else
                        new_map[i, j] = sum/count
                    end
                    
                end
            end
        end
    end
    return new_map
end