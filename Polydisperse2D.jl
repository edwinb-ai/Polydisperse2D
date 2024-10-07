using Random
using StaticArrays
using LinearAlgebra: dot
using DelimitedFiles: writedlm, readdlm
using Printf
using FastPow
using Statistics: mean
using Distributions: Gamma, Normal
using CodecZlib
using CellListMap
import CellListMap: copy_output, reset_output!, reducer
using Packmol: pack_monoatomic!

include("initialization.jl")
include("potentials.jl")
include("thermostat.jl")
include("pairwise.jl")
include("io.jl")

# Some numerical constants
const b_param = 1.0204081632653061
const a_param = 134.5526623421209

struct Parameters
    ρ::Float64
    ktemp::Float64
    n_particles::Int
end

function integrate_half(positions, velocities, forces, dt, boxl; pbc=true)
    # ! Important: There is a mass in the force term
    new_velocities = @. velocities + (forces * dt / 2.0)
    new_positions = @. positions + (new_velocities * dt)
    # Periodic boundary conditions
    if pbc
        new_positions = @. new_positions - boxl * round(new_positions / boxl)
    end

    return new_positions, new_velocities
end

function simulation(params::Parameters, pathname; eq_steps=100, prod_steps=500)
    rng = Random.Xoshiro()

    # For a polydisperse mixture, we need the diameters before anything else
    polydispersity = 0.11
    diameters = initialize_diameters(params.n_particles, polydispersity)

    # Now we compute the effective size of the box
    boxl = sqrt(sum(diameters .^ 2) / params.ρ)
    volume = boxl^2
    cutoff = 1.1
    # Set the timestep and the damping constant for the thermostat
    dt = 0.0001
    τ = 100.0 * dt
    # The degrees of freedom
    # Spatial dimension, in this case 2D simulations
    dimension = 2.0
    nf = dimension * (params.n_particles - 1.0)

    # Variables to accumulate results
    virial = 0.0
    nprom = 0
    kinetic_energy = 0.0

    # Initialize the system in a lattice configuration
    (system, diameters) = init_system(
        boxl, cutoff, pathname, diameters; n_particles=params.n_particles
    )
    # Initialize the velocities of the system by having the correct temperature
    velocities = initialize_velocities(params.ktemp, nf, rng, params.n_particles)
    # Adjust the particles using the velocities
    for i in eachindex(system.positions)
        system.positions[i] = @. system.positions[i] - (velocities[i] * dt)
    end
    # Zero out the arrays
    reset_output!(system.energy_and_forces)

    # Open files for trajectory and other things
    trajectory_file = joinpath(pathname, "production.xyz")
    eq_trajectory_file = joinpath(pathname, "equilibration.xyz")
    thermo_file = joinpath(pathname, "thermo.txt")
    eq_thermo_file = joinpath(pathname, "eq_thermo.txt")

    for step in 1:(eq_steps + prod_steps)
        for i in eachindex(system.positions, system.energy_and_forces.forces, velocities)
            f = system.energy_and_forces.forces[i]
            x = system.positions[i]
            v = velocities[i]
            (new_x, new_v) = integrate_half(x, v, f, dt, boxl)
            velocities[i] = new_v
            system.positions[i] = new_x
        end

        # Zero out arrays
        reset_output!(system.energy_and_forces)
        # Compute energy and forces
        map_pairwise!(
            (x, y, i, j, d2, output) ->
                energy_and_forces!(x, y, i, j, d2, output, diameters),
            system,
        )

        # Second half of the integration
        for i in eachindex(velocities, system.energy_and_forces.forces)
            f = system.energy_and_forces.forces[i]
            v = velocities[i]
            velocities[i] = @. v + (f * dt / 2.0)
        end

        # Always apply the thermostat and compute the kinetic energy
        bussi!(velocities, params.ktemp, nf, dt, τ, rng)
        kinetic_energy = compute_kinetic(velocities)

        # Accumulate the values of the virial for computing the pressure
        if mod(step, 100) == 0
            virial += system.energy_and_forces.virial
            nprom += 1
        end

        # Every few steps we save thermodynamic quantities to disk
        if mod(step, 1_000) == 0 && step > eq_steps
            ener_part = system.energy_and_forces.energy
            ener_part /= params.n_particles
            temperature = 2.0 * kinetic_energy / nf
            pressure = virial / (dimension * nprom * volume)
            pressure += params.ρ * temperature
            open(thermo_file, "a") do io
                writedlm(io, [step ener_part temperature pressure], " ")
            end
        end

        # Every few steps we save thermodynamic quantities to disk, equilibration
        if mod(step, 1_000) == 0 && step <= eq_steps
            ener_part = system.energy_and_forces.energy
            ener_part /= params.n_particles
            temperature = 2.0 * kinetic_energy / nf
            pressure = virial / (dimension * nprom * volume)
            pressure += params.ρ * temperature
            open(eq_thermo_file, "a") do io
                writedlm(io, [step ener_part temperature pressure], " ")
            end
        end

        # Save to disk the positions during equilibration
        if mod(step, 10_000) == 0 && step <= eq_steps
            # Write to file
            write_to_file(
                eq_trajectory_file,
                step,
                boxl,
                params.n_particles,
                system.positions,
                diameters,
            )
        end

        # Save to disk the positions during production
        if mod(step, 1_000) == 0 && step > eq_steps
            # Write to file
            write_to_file(
                trajectory_file, step, boxl, params.n_particles, system.positions, diameters
            )
        end
    end

    # Write the final configuration
    final_configuration = joinpath(pathname, "final.xyz")
    write_to_file(
        final_configuration, step, boxl, params.n_particles, system.positions, diameters; mode="w"
    )

    # Close all opened files
    close(thermo_file)
    # Also compress the trajectory file
    compress_gz(trajectory_file)

    return nothing
end

function main()
    densities = [0.95]
    ktemp = 1.4671
    n_particles = 2^14

    for d in densities
        params = Parameters(d, ktemp, n_particles)
        # Create a new directory with these parameters
        pathname = joinpath(@__DIR__, "N=$(n_particles)_density=$(@sprintf("%.4g", d))")
        # mkpath(pathname)
        # simulation(params, pathname; eq_steps=10_000_000, prod_steps=100_000_000)
        read_file(joinpath(pathname, "final.xyz"))
    end

    return nothing
end

main()
