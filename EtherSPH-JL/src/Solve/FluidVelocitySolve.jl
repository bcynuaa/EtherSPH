#=
  @ author: bcynuaa
  @ date: 2024-01-03 16:07:45
  @ description:
 =#

 function eachStep!(
    fluid_particles::FluidParticlesType where FluidParticlesType <: AbstractVector{<:FluidParticle},
    velocity_particles::WallParticlesType where WallParticlesType <: AbstractVector{<:VelocityParticle},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel,
    dr_forward_euler::DensityReinitializedForwardEuler,
    step::IntType where IntType <: Integer,
    check::Function
)::Nothing
    f_neighbours = findNeighbours(fluid_particles, smooth_kernel);
    f_w_neighbours = findNeighbours(fluid_particles, velocity_particles, smooth_kernel);
    Threads.@threads for neighbour in f_neighbours
        continuity!(fluid_particles[neighbour.i_], fluid_particles[neighbour.j_], neighbour);
    end
    Threads.@threads for i_f_particle in eachindex(fluid_particles)
        updateDensity!(fluid_particles[i_f_particle], dr_forward_euler.dt_);
        updatePressure!(fluid_particles[i_f_particle], wc_lm);
    end
    Threads.@threads for neighbour in f_neighbours
        pressureForce!(fluid_particles[neighbour.i_], fluid_particles[neighbour.j_], neighbour, smooth_kernel);
        viscosityForce!(fluid_particles[neighbour.i_], fluid_particles[neighbour.j_], neighbour, smooth_kernel, wc_lm);
    end
    Threads.@threads for neighbour in f_w_neighbours
        wallForce!(fluid_particles[neighbour.i_], velocity_particles[neighbour.j_], neighbour, smooth_kernel);
        pressureForce!(fluid_particles[neighbour.i_], neighbour, smooth_kernel, wc_lm);
        # viscosityForce!(fluid_particles[neighbour.i_], neighbour, smooth_kernel, wc_lm);
        viscosityForce!(fluid_particles[neighbour.i_], velocity_particles[neighbour.j_], neighbour, smooth_kernel, wc_lm);
    end
    Threads.@threads for i_f_particle in eachindex(fluid_particles)
        updateVelocity!(fluid_particles[i_f_particle], dr_forward_euler.dt_, wc_lm.body_force_vec_);
        updatePosition!(fluid_particles[i_f_particle], dr_forward_euler.dt_);
    end
    if isDensityReinitializedStep(step, dr_forward_euler)
        reconstructScalar!(fluid_particles, :rho_, f_neighbours, smooth_kernel);
        out_of_bounds_index = [];
        for i_f_particle in eachindex(fluid_particles)
            if !check(fluid_particles[i_f_particle])
                push!(out_of_bounds_index, i_f_particle);
            end
        end
        deleteat!(fluid_particles, out_of_bounds_index);
    end
    return nothing;
end

function solve!(
    fluid_particles::FluidParticlesType where FluidParticlesType <: AbstractVector{<:FluidParticle},
    velocity_particles::VelocityParticlesType where VelocityParticlesType <: AbstractVector{<:VelocityParticle},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel,
    dr_forward_euler::DensityReinitializedForwardEuler,
    vtp_io::VTPIO,
    check::Function
)::Nothing
    assureDirPathExist(vtp_io);
    writeVTP(0, 0., vtp_io, [fluid_particles, velocity_particles]);
    for step in ProgressBar(1: dr_forward_euler.total_step_)
        # eachStep!(fluid_particles, velocity_particles, smooth_kernel, wc_lm, dr_forward_euler, step, check);
        try
            eachStep!(fluid_particles, velocity_particles, smooth_kernel, wc_lm, dr_forward_euler, step, check);
        catch e
            println(e);
            println("step: ", step);
            writeVTP(step, step * dr_forward_euler.dt_, vtp_io, [fluid_particles, velocity_particles]);
            break;
        end
        if isOutputStep(step, dr_forward_euler)
            writeVTP(div(step, dr_forward_euler.output_step_), step * dr_forward_euler.dt_, vtp_io, [fluid_particles, velocity_particles]);
        end
    end
    return nothing;
end