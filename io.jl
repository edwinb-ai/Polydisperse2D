function write_to_file(filepath, step, boxl, n_particles, positions, diameters)
    # Write to file
    open(filepath, "a") do io
        println(io, n_particles)
        Printf.@printf(
            io,
            "Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 0.0\" Properties=type:I:1:id:I:1:radius:R:1:pos:R:2 Time=%d\n",
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