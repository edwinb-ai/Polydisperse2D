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
using Packmol: pack_monoatomic!

# Some numerical constants
const b_param = 1.0204081632653061
const a_param = 134.5526623421209

struct Parameters
    ρ::Float64
    ktemp::Float64
    n_particles::Int
    polydispersity::Float64
end

include("initialization.jl")
include("potentials.jl")
include("thermostat.jl")
include("pairwise.jl")
include("io.jl")
include("minimize.jl")

function integrate_half!(positions, velocities, forces, dt::Float64, boxl::Float64)
    for i in eachindex(positions, forces, velocities)
        f = forces[i]
        x = positions[i]
        v = velocities[i]
        # ! Important: There is a mass in the force term
        velocities[i] = @. v + (f * dt / 2.0)
        positions[i] = @. x + (velocities[i] * dt)
        positions[i] = @. positions[i] - boxl * round(positions[i] / boxl)
    end

    return nothing
end

function integrate_second_half!(velocities, forces, dt)
    for i in eachindex(velocities, forces)
        f = forces[i]
        v = velocities[i]
        velocities[i] = @. v + (f * dt / 2.0)
    end

    return nothing
end

function simulation(
    params::Parameters, pathname; from_file="", eq_steps=100, prod_steps=500
)
    rng = Random.Xoshiro()

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
    kinetic_temperature = 0.0

    # Initialize the system
    (system, diameters, volume, boxl) = initialize_simulation(
        params, pathname; polydispersity=params.polydispersity, file=from_file
    )

    # Initialize the velocities of the system by having the correct temperature
    velocities = initialize_velocities(params.ktemp, nf, rng, params.n_particles)
    # Adjust the particles using the velocities
    for i in eachindex(system.positions)
        system.positions[i] = @. system.positions[i] - (velocities[i] * dt)
    end
    # Zero out the arrays
    reset_output!(system.energy_and_forces)

    # Remove the files if they existed, and return the files handles
    (trajectory_file, eq_trajectory_file, thermo_file, eq_thermo_file) = open_files(
        pathname
    )
    # Formatting string
    format_string = Printf.Format("%d %.6f %.6f %.6f\n")

    # The main loop of the simulation
    for step in 1:(eq_steps + prod_steps)
        # First half of the integration
        integrate_half!(
            system.positions, velocities, system.energy_and_forces.forces, dt, boxl
        )

        # Zero out arrays
        reset_output!(system.energy_and_forces)
        # Compute energy and forces
        CellListMap.map_pairwise!(
            (x, y, i, j, d2, output) ->
                energy_and_forces!(x, y, i, j, d2, output, diameters),
            system,
        )

        # Second half of the integration
        integrate_second_half!(velocities, system.energy_and_forces.forces, dt)

        # Always apply the thermostat and compute the kinetic energy
        bussi!(velocities, params.ktemp, nf, dt, τ, rng)
        kinetic_energy = compute_kinetic(velocities)
        temperature = 2.0 * kinetic_energy / nf

        # Accumulate the values of the virial for computing the pressure
        if mod(step, 10) == 0
            virial += system.energy_and_forces.virial
            kinetic_temperature += temperature
            nprom += 1
        end

        # Every few steps we save thermodynamic quantities to disk
        if mod(step, 1_000) == 0 && step > eq_steps
            ener_part = system.energy_and_forces.energy
            ener_part /= params.n_particles
            average_temperature = kinetic_temperature / nprom
            pressure = virial / (dimension * nprom * volume)
            pressure += params.ρ * average_temperature
            open(thermo_file, "a") do io
                Printf.format(
                    io, format_string, step, ener_part, average_temperature, pressure
                )
            end

            # Reset the values of accumulators
            virial = 0.0
            kinetic_temperature = 0.0
            nprom = 0.0
        end

        # Every few steps we save thermodynamic quantities to disk, equilibration
        if mod(step, 10_000) == 0 && step <= eq_steps
            ener_part = system.energy_and_forces.energy
            ener_part /= params.n_particles
            average_temperature = kinetic_temperature / nprom
            pressure = virial / (dimension * nprom * volume)
            pressure += params.ρ * average_temperature
            open(eq_thermo_file, "a") do io
                Printf.format(
                    io, format_string, step, ener_part, average_temperature, pressure
                )
            end

            # Reset the values of accumulators
            virial = 0.0
            kinetic_temperature = 0.0
            nprom = 0.0
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
    step = eq_steps + prod_steps
    write_to_file(
        final_configuration,
        step,
        boxl,
        params.n_particles,
        system.positions,
        diameters;
        mode="w",
    )

    # Also do a minimization of the system
    minimize(system, diameters, boxl; tolerance=1e-5)
    min_configuration = joinpath(pathname, "minimized.xyz")
    write_to_file(
        min_configuration,
        step,
        boxl,
        params.n_particles,
        system.positions,
        diameters;
        mode="w",
    )

    # Also compress the trajectory file, if it exists
    if isfile(trajectory_file)
        compress_gz(trajectory_file)
    end

    if isfile(eq_trajectory_file)
        compress_gz(eq_trajectory_file)
    end

    return nothing
end

function read_minimize(params::Parameters, pathname; from_file="")
    # Initialize the system
    (system, diameters, _, boxl) = initialize_simulation(
        params, pathname; polydispersity=params.polydispersity, file=from_file
    )

    # Zero out arrays
    reset_output!(system.energy_and_forces)
    # Compute energy and forces
    CellListMap.map_pairwise!(
        (x, y, i, j, d2, output) -> energy_and_forces!(x, y, i, j, d2, output, diameters),
        system,
    )

    current_energy = system.energy_and_forces.energy
    println("Initial energy: $(current_energy)")

    minimize(system, diameters, boxl)

    # Zero out arrays
    reset_output!(system.energy_and_forces)
    # Compute energy and forces
    CellListMap.map_pairwise!(
        (x, y, i, j, d2, output) -> energy_and_forces!(x, y, i, j, d2, output, diameters),
        system,
    )

    current_energy = system.energy_and_forces.energy
    println("Final energy: $(current_energy)")

    final_configuration = joinpath(pathname, "minimized.xyz")
    step = 0
    write_to_file(
        final_configuration,
        step,
        boxl,
        params.n_particles,
        system.positions,
        diameters;
        mode="w",
    )

    return nothing
end

function main()
    densities = [0.955]
    ktemp = 1.4671
    n_particles = 2^12
    polydispersity = 0.15

    for d in densities
        params = Parameters(densities[1], ktemp, n_particles, polydispersity)
        # Create a new directory with these parameters
        pathname = joinpath(
            @__DIR__,
            "N=$(n_particles)_density=$(@sprintf("%.4g", d))_Δ=$(@sprintf("%.2g", polydispersity))",
        )
        mkpath(pathname)
        simulation(params, pathname; eq_steps=100_000, prod_steps=1)
        # read_minimize(params, pathname; from_file=joinpath(pathname, "final.xyz"))
    end

    return nothing
end

main()
