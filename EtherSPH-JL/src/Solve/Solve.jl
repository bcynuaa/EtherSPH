#=
  @ author: bcynuaa
  @ date: 2023-11-28 20:11:59
  @ description:
 =#

function eachStep!(
    particle_pool::ParticlePool,
    body_force_vec::AbstractVector,
    kernel::AbstractKernel,
    wc_liquid_model::WeakCompressibleLiquidModel,
    forward_euler::ForwardEuler
)::Nothing
    particle_neighbours = findNeighbours(particle_pool, kernel);
    # balance mass
    for neighbour in particle_neighbours
        balanceMass!(particle_pool.particles_[neighbour.i_], particle_pool.particles_[neighbour.j_], neighbour, wc_liquid_model);
    end
    # update density and pressure
    for particle in particle_pool.particles_
        updateDensity!(particle, forward_euler.dt_);
        updatePressure!(particle, wc_liquid_model);
    end
    # pressure force and viscosity force
    for neighbour in particle_neighbours
        pressureForce!(particle_pool.particles_[neighbour.i_], particle_pool.particles_[neighbour.j_], neighbour, wc_liquid_model);
        viscosityForce!(particle_pool.particles_[neighbour.i_], particle_pool.particles_[neighbour.j_], neighbour, kernel, wc_liquid_model);
    end
    # update velocity and position
    for particle in particle_pool.particles_
        updateVelocity!(particle, forward_euler.dt_, body_force_vec);
        updatePosition!(particle, forward_euler.dt_);
    end
    return nothing;
end

function solve!(
    particle_pool::ParticlePool,
    body_force_vec::AbstractVector,
    kernel::AbstractKernel,
    wc_liquid_model::WeakCompressibleLiquidModel,
    forward_euler::ForwardEuler,
    vtp_io::VTPIO
)::Nothing
    # assure dir path exist
    assureDirPathExist(vtp_io);
    # write initial state
    writeVtp(particle_pool, 0, vtp_io);
    # start time loop
    for step in ProgressBar(1: forward_euler.total_step_)
        # each step
        eachStep!(particle_pool, body_force_vec, kernel, wc_liquid_model, forward_euler);
        # write vtp file
        if isStepForOutput(forward_euler, step)
            writeVtp(particle_pool, div(step, forward_euler.output_step_), vtp_io);
        end
    end
    return nothing;
end