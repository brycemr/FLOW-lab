using LinearAlgebra

gamma = [0; 0; 1]
d = 2.0
timestep = 0.01
p1 = [0.0, -d/2, 0.0]
p2 = [0.0, d/2, 0.0]
p3 = [1.0, d/2, 0.0]
p4 = [1.0, -d/2, 0.0]
currentVelocities = zeros(4, 3)

"""
velocities(gamma, P1, P2, P3, P4)

calculates the resultant velocities for each vortex given the current vortex positions and value for gamma
"""
function velocities(gamma, P1, P2, P3, P4)
    r12 = P1 - P2
    r13 = P1 - P3
    r14 = P1 - P4
    r23 = P2 - P3
    r24 = P2 - P4
    r34 = P3 - P4
    v1 = (cross(gamma,r12) / (2*pi*(norm(r12)^2))) + (cross(gamma,r13) / (2*pi*(norm(r13)^2))) + (cross((gamma.*-1),r14) / (2*pi*(norm(r14)^2)))
    v2 = (cross((gamma.*-1),(r12.*-1)) / (2*pi*(norm(r12)^2))) + (cross(gamma,r23) / (2*pi*(norm(r23)^2))) + (cross((gamma.*-1),r24) / (2*pi*(norm(r24)^2)))
    v3 = (cross((gamma.*-1),(r13.*-1)) / (2*pi*(norm(r13)^2))) + (cross(gamma,(r23.*-1)) / (2*pi*(norm(r23)^2))) + (cross((gamma.*-1),r34) / (2*pi*(norm(r34)^2)))
    v4 = (cross((gamma.*-1),(r14.*-1)) / (2*pi*(norm(r14)^2))) + (cross(gamma,(r24.*-1)) / (2*pi*(norm(r24)^2))) + (cross(gamma,(r34.*-1)) / (2*pi*(norm(r34)^2)))
    velocities = hcat(v1, v2, v3, v4)
    return velocities
end

using Plots
allp1 = p1
allp2 = p2
allp3 = p3
allp4 = p4
for i in 0:timestep:40
    global currentVelocities = velocities(gamma, p1, p2, p3, p4)
    #println(currentVelocities[:,1])
    global p1 = p1 + currentVelocities[:,1].*timestep
    global p2 = p2 + currentVelocities[:,2].*timestep
    global p3 = p3 + currentVelocities[:,3].*timestep
    global p4 = p4 + currentVelocities[:,4].*timestep
    global allp1 = hcat(allp1,p1)
    global allp2 = hcat(allp2,p2)
    global allp3 = hcat(allp3,p3)
    global allp4 = hcat(allp4,p4)
end

plot(allp1[1,:],allp1[2,:],
    size = (1800, 600),
    lw = 2, lc = :blue,
    grid = false,
    legend = false,
    yticks = [0])

plot!(allp2[1,:],allp2[2,:], lw = 2, lc = :blue)
plot!(allp3[1,:],allp3[2,:], lw = 2, lc = :red)
plot!(allp4[1,:],allp4[2,:], lw = 2, lc = :red)
