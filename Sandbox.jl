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