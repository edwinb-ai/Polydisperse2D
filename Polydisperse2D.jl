using Random
using StaticArrays
using LinearAlgebra
using DelimitedFiles: writedlm, readdlm
using Printf
using FastPow
using Statistics: mean
using Distributions: Gamma, Normal
using CodecZlib
using CellListMap
import CellListMap: copy_output, reset_output!, reducer

struct Parameters
    ρ::Float64
    n_particles::Int
    polydispersity::Float64
    dt::Float64
end

include("initialization.jl")
include("potentials.jl")
include("thermostat.jl")
include("pairwise.jl")
include("io.jl")
include("minimize.jl")
include("integrate.jl")
include("simulation.jl")

function main()
    # Read the density from the command line
    density = 0.99
    ktemp = 1.4671
    n_particles = 2^10
    polydispersity = 0.15
    dt = 0.0001

    params = Parameters(density, n_particles, polydispersity, dt)
    pathname = joinpath(
        @__DIR__,
        "test_N=$(n_particles)_density=$(@sprintf("%.4g", density))_Δ=$(@sprintf("%.2g", polydispersity))",
    )
    mkpath(pathname)

    # We initialize the state of the simulation
    thermostat = NVT(ktemp, 100.0 * dt)
    state = initialize_state(params, ktemp, pathname)
    run_simulation!(state, params, thermostat, 1_000_000, 100_000, pathname; compress=true)

    # Now we do NVE
    run_simulation!(
        state,
        params,
        NVE(),
        1_000_000,
        100_000,
        pathname;
        traj_name="production.xyz",
        thermo_name="production_thermo.txt",
        log_times=true,
        compress=true,
    )

    return nothing
end

main()
