function minimize(system, diameters, boxl; dt=0.0001, tolerance=1e-4)
    # FIRE Parameters
    α_start = 0.1
    α = α_start
    N_min = 5
    f_inc = 1.1
    f_dec = 0.5
    α_dec = 0.99
    dt_max = 10.0 * dt
    N_positive = 0

    n_particles = length(diameters)
    velocities = [StaticArrays.@SVector zeros(2) for _ in 1:n_particles]

    converged = false
    while !converged
        # Zero out arrays
        reset_output!(system.energy_and_forces)
        # Compute energy and forces
        CellListMap.map_pairwise!(
            (x, y, i, j, d2, output) ->
                energy_and_forces!(x, y, i, j, d2, output, diameters),
            system,
        )
        # Update velocities
        for i in eachindex(system.energy_and_forces.forces)
            f = system.energy_and_forces.forces[i]
            v = velocities[i]
            velocities[i] = @. v + f * dt
        end

        # Compute the power
        power = 0.0
        for i in eachindex(velocities)
            f = system.energy_and_forces.forces[i]
            v = velocities[i]
            power += sum(dot(f, v))
        end

        if power > 0.0
            for i in eachindex(system.energy_and_forces.forces)
                f = system.energy_and_forces.forces[i]
                v = velocities[i]
                v_norm = norm(v)
                f_norm = norm(f)
                if f_norm != 0 && v_norm != 0
                    velocities[i] = @. (1.0 - α) * v + α * (v_norm / f_norm) * f
                end
            end
            N_positive += 1
            if N_positive > N_min
                dt = min(dt * f_inc, dt_max)
                α *= α_dec
            end
        else
            # Reset velocities and parameters
            for i in eachindex(velocities)
                velocities[i] = StaticArrays.@SVector zeros(2)
            end
            dt *= f_dec
            α = α_start
            N_positive = 0
        end

        # Update positions
        for i in eachindex(system.positions)
            r = system.positions[i]
            v = velocities[i]
            system.positions[i] = @. r + v * dt
            # Apply PBC
            system.positions[i] = @. system.positions[i] -
                boxl * round(system.positions[i] / boxl)
        end

        # Convergence Check
        max_force = maximum(norm.(system.energy_and_forces.forces))
        if max_force < tolerance
            converged = true
            if converged
                @info "Converged."
                @info max_force
            end
        end
    end

    return nothing
end