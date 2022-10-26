input = [1.0, 7.0, 5.0, 4.0, 3.0, 1.0]
nturbines = 4

adjacency = zeros((nturbines,nturbines))
println(adjacency)

#lengths = ff.turbine_spacing(turbinex,turbiney)
global shift = 2
for i = 1:nturbines
    for j = i+1:nturbines
        #The adjacency matrix is symmetric
        adjacency[i,j] = input[i+j-shift]
        adjacency[j,i] = input[i+j-shift]
    end
    global shift = 1
end
println(adjacency)
total = 0
visitedTurbines = [1]
global nextTurbine = 1
while length(visitedTurbines) < nturbines
    global nextTurbine += 1
    min = Inf
    println(visitedTurbines, nextTurbine)
    for turbine in visitedTurbines
        if adjacency[turbine,nextTurbine] < min
            min = adjacency[turbine,nextTurbine]
        end
    end
    push!(visitedTurbines, nextTurbine)
    println(min)
    global total += min
end     

println(total)