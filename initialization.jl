function initialize_lattice(n_particles::Int, boxl::Float64)
    # Determine the size of the lattice (side length of the square)
    side_length = ceil(Int, sqrt(n_particles))  # Round up to ensure all particles fit

    # Calculate the spacing based on the box size
    spacing = boxl / side_length  # Ensure particles fit in the box

    # Initialize arrays to hold the x and y positions of the particles
    x_positions = []
    y_positions = []

    # Create the lattice by filling in the grid
    for i in 1:side_length
        for j in 1:side_length
            if length(x_positions) < n_particles  # Only place up to n_particles particles
                # Scale positions based on the spacing and shift to center particles in the box
                push!(x_positions, (i - 0.5) * spacing)
                push!(y_positions, (j - 0.5) * spacing)
            else
                break
            end
        end
    end

    # Return the particle positions as a tuple of arrays
    return x_positions, y_positions
end

function initialize_random(unitcell, npart, rng; tol=1.0)
    coordinates = unitcell[1] * rand(rng, StaticArrays.SVector{2,Float64}, npart)
    pack_monoatomic!(coordinates, unitcell, tol; parallel=false, iprint=100)

    return coordinates
end

function initialize_diameters(n_particles, polydispersity)
    # Now generate an array of normally distributed diameters
    mean_diameter = 1.0
    diameters = zeros(n_particles)
    poly_dist = Normal(mean_diameter, polydispersity)

    for i in 1:n_particles
        new_diameter = rand(poly_dist)
        # Make sure no diameters are larger than 1
        while new_diameter > 1.0
            new_diameter = rand(poly_dist)
        end
        diameters[i] = new_diameter
    end

    @info "Mean diameter: $(mean(diameters))"

    return diameters
end

function init_system(boxl, cutoff, pathname, diameters; n_particles=2^8)
    # Define the size of the simulation box for the cell lists
    unitcell = [boxl, boxl]

    # We can create normal arrays for holding the particles' positions
    (x, y) = initialize_lattice(n_particles, boxl)

    # Now we change the arrays to static versions of it
    positions = [@SVector [i, j] for (i, j) in zip(x, y)]

    # Save the initial configuration to a file
    filepath = joinpath(pathname, "initial.xyz")
    write_to_file(filepath, 0, boxl, n_particles, positions, diameters)

    # Initialize system
    system = CellListMap.ParticleSystem(;
        xpositions=positions,
        unitcell=[boxl, boxl],
        cutoff=cutoff,
        output=EnergyAndForces(0.0, 0.0, similar(positions)),
        output_name=:energy_and_forces,
        parallel=false,
    )

    return system, diameters
end

function initialize_simulation(params::Parameters, pathname; polydispersity=0.11, file="")
    # Always leave a fixed cutoff
    cutoff = 1.1
    diameters = []
    positions = []
    boxl = 0.0
    volume = 0.0
    system = nothing

    if isfile(file)
        @info "Reading from file..."
        (boxl, positions, diameters) = read_file(file)
        volume = boxl^2

        # Initialize system
        system = CellListMap.ParticleSystem(;
            xpositions=positions,
            unitcell=[boxl, boxl],
            cutoff=cutoff,
            output=EnergyAndForces(0.0, 0.0, similar(positions)),
            output_name=:energy_and_forces,
            parallel=false,
        )

        # Save the initial configuration to a file
        filepath = joinpath(pathname, "initial.xyz")
        write_to_file(filepath, 0, boxl, params.n_particles, positions, diameters)

        @info "Reading done."
    else
        @info "Creating a new system."
        # For a polydisperse mixture, we need the diameters before anything else
        diameters = initialize_diameters(params.n_particles, polydispersity)

        # Now we compute the effective size of the box
        boxl = sqrt(sum(diameters .^ 2) / params.ρ)
        volume = boxl^2

        # Initialize the system in a lattice configuration
        (system, diameters) = init_system(
            boxl, cutoff, pathname, diameters; n_particles=params.n_particles
        )
    end

    return system, diameters, volume, boxl
end

function initialize_velocities(ktemp, nf, rng, n_particles)
    # Initilize the random numbers of the velocities
    velocities = [StaticArrays.@SVector zeros(2) for _ in 1:n_particles]
    sum_v = StaticArrays.@MVector zeros(2)
    sum_v2 = 0.0

    for i in eachindex(velocities)
        velocities[i] = randn(rng, size(velocities[i]))
        # Collect the center of mass, momentum = 1
        sum_v .+= velocities[i]
    end

    sum_v ./= n_particles

    for i in eachindex(velocities)
        # Remove the center of mass momentum
        velocities[i] = velocities[i] .- sum_v
        sum_v2 += sum(abs2, velocities[i])
    end

    fs = sqrt(ktemp / (sum_v2 / nf))
    for i in eachindex(velocities)
        velocities[i] = velocities[i] .* fs
    end

    return velocities
end