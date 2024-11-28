struct SimulationState{T,U,V}
    # This field contains the cell lists for the system itself
    system::T
    # The array that contains the diameters of the particles
    diameters::U
    # The RNG
    rng::V
    # The size of the simulation box
    boxl::Float64
    # The container for the velocities
    velocities::Vector{SVector{2,Float64}}
    # The images for the particles
    images::Vector{SVector{2,Float64}}
end

function initialize_state(params::Parameters, ktemp::Float64, pathname; from_file="")
    rng = Random.Xoshiro()

    # The degrees of freedom
    # Spatial dimension, in this case 2D simulations
    dimension = 2.0
    nf = dimension * (params.n_particles - 1.0)

    # Initialize the system
    (system, diameters, boxl) = initialize_simulation(
        params, pathname; polydispersity=params.polydispersity, file=from_file
    )
    # Initialize the velocities of the system by having the correct temperature
    velocities = initialize_velocities(ktemp, nf, rng, params.n_particles)
    # Adjust the particles using the velocities
    for i in eachindex(system.positions)
        system.positions[i] = @. system.positions[i] - (velocities[i] * params.dt)
    end
    # Zero out the arrays
    reset_output!(system.energy_and_forces)

    # Initialize the array of images
    images = [@SVector zeros(dimension) for _ in eachindex(velocities)]

    state = SimulationState(system, diameters, rng, boxl, velocities, images)

    return state
end

function ensemble_step!(
    ::NVE, velocities, params::Parameters, nf::Float64, state::SimulationState
)
    kinetic_energy = compute_kinetic(velocities)
    temperature = 2.0 * kinetic_energy / nf
    return temperature
end

function ensemble_step!(
    ensemble::NVT, velocities, params::Parameters, nf::Float64, state::SimulationState
)
    # Apply thermostat, e.g., Bussi thermostat
    bussi!(velocities, ensemble.ktemp, nf, params.dt, ensemble.tau, state.rng)

    kinetic_energy = compute_kinetic(velocities)
    temperature = 2.0 * kinetic_energy / nf
    return temperature
end

function finalize_simulation!(
    trajectory_file::String,
    pathname::String,
    total_steps::Int,
    state::SimulationState,
    params::Parameters,
)
    final_configuration = joinpath(pathname, "final.xyz")
    write_to_file(
        final_configuration,
        total_steps,
        state.boxl,
        params.n_particles,
        state.system.positions,
        state.diameters;
        mode="w",
    )

    if isfile(trajectory_file)
        compress_gz(trajectory_file)
    end
end

function run_simulation!(
    state::SimulationState,
    params::Parameters,
    ensemble::Ensemble,
    total_steps::Int,
    frequency::Int,
    pathname;
    traj_name="trajectory.xyz",
    thermo_name="thermo.txt",
)
    # Remove the files if they existed, and return the files handles
    (trajectory_file, thermo_file) = open_files(pathname, traj_name, thermo_name)
    format_string = Printf.Format("%d %.6f %.6f %.6f\n")

    # Extract parameters from the state
    system = state.system
    velocities = state.velocities
    diameters = state.diameters
    pbc = true

    # Degrees of freedom
    dimension = 2.0
    nf = dimension * (params.n_particles - 1.0)

    # Compute the volume
    volume = state.boxl^dimension

    # Variables to accumulate results
    virial = 0.0
    nprom = 0
    kinetic_temperature = 0.0

    for step in 1:total_steps
        # Perform integration
        integrate_half!(
            system.positions,
            velocities,
            system.energy_and_forces.forces,
            params.dt,
            state.boxl;
            pbc=pbc,
        )
        reset_output!(system.energy_and_forces)
        CellListMap.map_pairwise!(
            (x, y, i, j, d2, output) ->
                energy_and_forces!(x, y, i, j, d2, output, diameters),
            system,
        )
        integrate_second_half!(velocities, system.energy_and_forces.forces, params.dt)

        # Apply ensemble-specific logic
        temperature = ensemble_step!(ensemble, velocities, params, nf, state)

        # Accumulate values for thermodynamics
        if mod(step, 10) == 0
            virial += system.energy_and_forces.virial
            kinetic_temperature += temperature
            nprom += 1
        end

        # Output thermodynamic quantities periodically
        if mod(step, frequency) == 0
            ener_part = system.energy_and_forces.energy / params.n_particles
            avg_temp = kinetic_temperature / nprom
            pressure = virial / (dimension * nprom * volume) + params.ρ * avg_temp
            open(thermo_file, "a") do io
                Printf.format(io, format_string, step, ener_part, avg_temp, pressure)
            end
            virial, kinetic_temperature, nprom = 0.0, 0.0, 0
        end

        # Output trajectory periodically
        if mod(step, frequency) == 0
            write_to_file(
                trajectory_file,
                step,
                state.boxl,
                params.n_particles,
                system.positions,
                diameters,
            )
        end
    end

    # Final output and cleanup
    finalize_simulation!(trajectory_file, pathname, total_steps, state, params)
    return nothing
end