# Learning Julia
# FLOW Lab LTRAD Program
# Author : Bryce Richard
# Date : 5/19/2022
# Introductory program for learning basics of Julia

using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Plots")
using CSV
using DataFrames
using Plots

a = 2 + 2
b = 3
println("The coolest thing is that $b is not equal to $a")

"""
nacathickness(x, m)

nacathickness takes in a maximum thickness (m) and x-position (x) and returns the airfoil thickness at that position

"""
function nacathickness(x, m)
    output = 10*m*(0.2969*sqrt(x) - 0.126*x - 0.3537*(x^2)+0.2843*(x^3)-0.1015*(x^4))
    return output
end

#=
FOR LOOP to validate nacathickness

for i in 0:10
    x = round(i*0.10, digits=1)
    t = nacathickness(x, 0.1)
    println("$x -> $t")
end
=#

"""
nacacamber(x, p, c)

nacathickness calculates the camber given an x position, the maximum camber of the airfoil, and the position of maximum camber

"""
function nacacamber(x, p, c)
    camber = (c*(2*p*x - (x^2)))/(p^2)
    if x == 0 && p == 0
        camber = 0
    end
    if x > p
        camber = (c*(1-2*p + 2*p*x - (x^2)))/((1 - p)^2)
    end
    return camber
end

#=
FOR LOOP to validate nacacamber

for i in 0:10
    x = round(i*0.10, digits=1)
    z = nacacamber(x, 0.4, 0.02)
    println("$x -> $z")
end
=#

x = []
precision = 0.01
for i in 0:precision:1
    append!(x,i)
end

println(x)
"""
nacacoordinates(c, p, t, x)

calculates the upper and lower coordinates along an airfoil given max camber (c), max camber postion (p), max thickness (t), and x positions

"""
function nacacoordinates(c::Int, p::Int, t::Int, x)
    halfthickness = broadcast(/, nacathickness.(x,t/100), 2)
    upper =  broadcast(+, nacacamber.(x,p/10,c/100) , halfthickness)
    lower =  broadcast(-, nacacamber.(x,p/10,c/100) , halfthickness)
    return vcat(reverse(lower),upper)
end

naca2412 = nacacoordinates(2, 4, 12, x)
naca0008 = nacacoordinates(0, 0, 8, x)
naca0012 = nacacoordinates(0, 0, 12, x)
naca4412 = nacacoordinates(4, 4, 12, x)
naca2424 = nacacoordinates(2, 4, 24, x)

xall = vcat(reverse(x),x)

#=
FOR LOOP to validate nacacoordinates and print out results

for i in 1:11
    currx = x[i]
    currupper = upper[i]
    currlower = lower[i]
    println("$currx | $currupper \t| $currlower")
end
=#

#=
HOW TO READ IN CSV FILE
df = CSV.read("E203.csv", DataFrame)
print(df)
=#

plot(xall,naca0008,
    size = (800, 400),
    grid = false,
    legend_font_pointsize = 12,
    foreground_color_legend = nothing,
    ylims = (-0.15,0.3),
    xticks = (0:1:1),
    yticks = [-0.15,0,0.25],
    label = "NACA 0008", lw = 4 )
plot!(xall,naca0012, label = "NACA 0012", lw = 3, ls = :dot)
plot!(xall,naca2412, label = "NACA 2412", lw = 3, lc = :purple)
plot!(xall,naca4412, label = "NACA 4412", lw = 3, lc = :green, ls = :dash)
plot!(xall,naca2424, label = "NACA 2424", lw = 3, ls = :dashdotdot)

filename = "Naca_Airfoils_Plot.pdf"
savefig(filename)