module simpleMD
# Load Packages
using CSV
using DataFrames
using LinearAlgebra
using ProgressMeter

#using Plots
#using Random
wrappedDispacement(A,B, L)    = abs.(abs.(mag(A-B).-L/2).-L/2)
mag(x::Vector) = sqrt(x·x)

mutable struct Particle # Define Particle Type
    # Charge
    q::Float64  # TODO vary by species
    # Mass
    m::Float64
    # Kinematic Variables
    # TODO define types and conversions for fields
    pos::Vector{Float64}
    vel::Vector{Float64}
    acc::Vector{Float64}
    acc0::Vector{Float64}   
    # Energy
    #K::Float64 # Calc. frm (m,v),(funct <- particle) =>K
    U::Float64
    # Inner Constructor Methods
    # Convenience constructor
    Particle(pos::Vector) = new(1.0, 1.0, pos, zeros(3), zeros(3),zeros(3), 0.0)
    Particle(pos::Vector, vel::Vector) =
        new(1.0, 1.0, pos, vel, zeros(3),zeros(3), 0.0, 0.0)
    Particle(q, m, pos) = new(q, m, pos, zeros(3), zeros(3),zeros(3), 0.0)
    Particle(q, m, pos::Vector, vel::Vector) =
        new(q, m, pos, vel, zeros(3), 0.5 * m * mag(vel)^2)
end

function apot(Pi::Particle, Pj::Particle, L::Number)
    # Returns a tuple containing the force::vector and energy of interaction
    disp = Pi.pos - Pj.pos  # displacment (w/ PBC)
    wrappedDispacement    = abs.(abs.(disp.-L/2).-L/2);
    dist                  = mag(wrappedDispacement);
    a = (Pi.q * Pj.q) * disp / (Pi.m*dist^3)     # Return calculated force
    U = (Pi.q * Pj.q) / dist     # Return calculated potential
    a, U
end

function step!(apot, Ps::Vector{Particle}, L::Float64; Δt=1e-6)
    # apot must be a function returning a tuple with the acceleration vector
    # in the first position and the interaction energy in the second.
   for Pi in Ps
        # Reset Acceleration & Interaction Potential
        Pi.acc0 = Pi.acc
        Pi.acc = zeros(3)
        Pi.U   = 0
        for Pj in Ps
            if Pi ≠ Pj
                #println(fpot(Pi, Pj, L))
                Pi.acc, Pi.U = (Pi.acc, Pi.U) .+ apot(Pi, Pj, L)
            end
        end
    end

    for Pi in Ps
           # Update Velocity
           Pi.vel = Pi.vel + 0.5 * (Pi.acc + Pi.acc0) * Δt   # Uses avg of acc over timestep
           # Update Position
           Pi.pos = Pi.pos + Pi.vel * Δt + 0.5 * (Pi.acc) * (Δt)^2.0 # Verlet
           #Perform PBC wrapping
           Pi.pos = mod.(Pi.pos,L);
    end

end

N=100;
L=10.0;
Ntimes=200;

initpos=Positions = L * rand(Float64, (N, 3));

#data = CSV.File("savefile.xyz") |> Tables.matrix
#kindata = convert.(Float64, data[:, 2:end])
#initpos = kindata[:, 1:3]
#initvel = kindata[:,4:6]

particles = [Particle(convert(Vector, pos)) for pos in eachrow(initpos)]

p = Progress(Ntimes; showspeed=true)

@time for t in 1:Ntimes
    step!(apot, particles, L)
    next!(p)
end

end