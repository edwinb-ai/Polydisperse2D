function integrate_half!(
    positions, velocities, forces, dt::Float64, boxl::Float64; pbc=true
)
    for i in eachindex(positions, forces, velocities)
        f = forces[i]
        x = positions[i]
        v = velocities[i]
        # ! Important: There is a mass in the force term
        velocities[i] = @. v + (f * dt / 2.0)
        positions[i] = @. x + (velocities[i] * dt)
        if pbc
            positions[i] = @. positions[i] - boxl * round(positions[i] / boxl)
        end
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