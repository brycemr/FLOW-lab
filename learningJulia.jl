# Learning Julia
# FLOW Lab LTRAD Program
# Author : Bryce Richard
# Date : 5/19/2022
# Introductory program for learning basics of Julia
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


