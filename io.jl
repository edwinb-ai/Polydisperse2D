function write_to_file(filepath, step, boxl, n_particles, positions, diameters; mode="a")
    # Write to file
    open(filepath, mode) do io
        println(io, n_particles)
        Printf.@printf(
            io,
            "Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0\" Properties=type:I:1:id:I:1:radius:R:1:pos:R:2 Time=%.6g\n",
            boxl,
            boxl,
            step
        )
        for i in eachindex(diameters, positions)
            particle = positions[i]
            # velocity = velocities[i]
            Printf.@printf(
                io,
                "%d %d %lf %lf %lf\n",
                1,
                i,
                diameters[i] / 2.0,
                particle[1],
                particle[2],
                # velocity[1],
                # velocity[2],
            )
        end
    end

    return nothing
end

function write_to_file_lammps(
    filepath, step, boxl, n_particles, positions, diameters; mode="w"
)
    # Write to file
    open(filepath, mode) do io
        Printf.@printf(io, "ITEM: TIMESTEP\n%d\n", step)
        Printf.@printf(io, "ITEM: NUMBER OF ATOMS\n%d\n", n_particles)
        Printf.@printf(
            io,
            "ITEM: BOX BOUNDS pp pp pp\n0.0 %lf\n0.0 %lf\n0.0 0.0\nITEM: ATOMS id type xu yu zu\n",
            boxl,
            boxl
        )
        for i in eachindex(diameters, positions)
            particle = positions[i]
            # velocity = velocities[i]
            Printf.@printf(
                io,
                "%d %d %lf %lf %lf\n",
                1,
                i,
                diameters[i] / 2.0,
                particle[1],
                particle[2],
                # velocity[1],
                # velocity[2],
            )
        end
    end

    return nothing
end

function read_file(filepath)
    # Initialize the variables with a default value
    n_particles = 0
    box_l = 0.0
    positions = []
    radii = []

    open(filepath, "r") do io
        # First line is the number of particles
        line = readline(io)
        n_particles = parse(Int64, line)

        # The second line we can skip for now
        line = split(readline(io), " ")
        # The fifth position included the box size length information
        box_l = parse(Float64, line[5])

        # Now read each line and gather the information
        positions = [@SVector(zeros(2)) for _ in 1:n_particles]
        radii = zeros(n_particles)

        for i in 1:n_particles
            line = split(readline(io), " ")
            parsed_line = parse.(Float64, line)
            positions[i] = SVector{2}(parsed_line[4:end])
            radii[i] = parsed_line[3]
        end
    end

    # The information read are the radii, so convert to diameters
    diameters = radii .* 2.0

    return box_l, positions, diameters
end

function compress_gz(filepath)
    # Attach the suffix to the original file
    output_file = filepath * ".gz"

    open(filepath, "r") do infile
        # Open the output file for writing, with gzip compression
        open(GzipCompressorStream, output_file, "w") do outfile
            # Write the contents of the input file to the compressed output file
            write(outfile, read(infile))
        end
    end

    # To avoid having double the files, we delete the original one
    rm(filepath)

    return nothing
end

function open_files(pathname, traj_name, thermo_name)
    # Open files for trajectory and other things
    trajectory_file = joinpath(pathname, traj_name)
    thermo_file = joinpath(pathname, thermo_name)

    files = [trajectory_file, thermo_file]

    for file in files
        if isfile(file)
            rm(file)
        end
    end

    return (trajectory_file, thermo_file)
end
