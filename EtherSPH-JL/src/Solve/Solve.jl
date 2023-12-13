#=
  @ author: bcynuaa
  @ date: 2023-12-06 19:20:56
  @ description:
 =#

using ProgressBars;

function eachStep!(
    fluid_particles::FluidParticlesType where FluidParticlesType <: AbstractVector{<:FluidParticle},
    wall_particles::WallParticlesType where WallParticlesType <: AbstractVector{<:WallParticle},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    xsph_wc_lm::XSPHWeaklyCompressibleLiquidModel,
    dr_forward_euler::DensityReinitializedForwardEuler,
    step::IntType where IntType <: Integer
)::Nothing
    f_neighbours = findNeighbours(fluid_particles, smooth_kernel);
    f_w_neighbours = findNeighbours(fluid_particles, wall_particles, smooth_kernel);
    Threads.@threads for neighbour in f_neighbours
        continuity!(fluid_particles[neighbour.i_], fluid_particles[neighbour.j_], neighbour, xsph_wc_lm);
    end
    Threads.@threads for i_f_particle in eachindex(fluid_particles)
        updateDensity!(fluid_particles[i_f_particle], dr_forward_euler.dt_);
        updatePressure!(fluid_particles[i_f_particle], xsph_wc_lm);
    end
    Threads.@threads for neighbour in f_neighbours
        pressureForce!(fluid_particles[neighbour.i_], fluid_particles[neighbour.j_], neighbour);
        viscosityForce!(fluid_particles[neighbour.i_], fluid_particles[neighbour.j_], neighbour, smooth_kernel, xsph_wc_lm);
    end
    Threads.@threads for neighbour in f_w_neighbours
        wallForce!(fluid_particles[neighbour.i_], wall_particles[neighbour.j_], neighbour, smooth_kernel, xsph_wc_lm);
    end
    Threads.@threads for i_f_particle in eachindex(fluid_particles)
        updateVelocity!(fluid_particles[i_f_particle], dr_forward_euler.dt_, xsph_wc_lm.body_force_vec_);
        updatePosition!(fluid_particles[i_f_particle], dr_forward_euler.dt_);
    end
    Threads.@threads for neighbour in f_neighbours
        positionCorrection!(fluid_particles[neighbour.i_], fluid_particles[neighbour.j_], dr_forward_euler.dt_, neighbour, xsph_wc_lm);
    end
    # Threads.@threads for neighbour in f_w_neighbours
    #     positionCorrection!(fluid_particles[neighbour.i_], dr_forward_euler.dt_, neighbour, xsph_wc_lm);
    # end
    if isDensityReinitializedStep(step, dr_forward_euler)
        reconstructScalar!(fluid_particles, :rho_, f_neighbours, smooth_kernel);
    end
    f_neighbours = nothing;
    f_w_neighbours = nothing;
    return nothing;
end

function solve!(
    particles_list::ParticlesListType where ParticlesListType <: AbstractVector{<:AbstractVector},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    xsph_wc_lm::XSPHWeaklyCompressibleLiquidModel,
    dr_forward_euler::DensityReinitializedForwardEuler,
    vtp_io::VTPIO
)::Nothing
    assureDirPathExist(vtp_io);
    writeVTP(0, 0., vtp_io, particles_list);
    for step in ProgressBar(1: dr_forward_euler.total_step_)
        try
            eachStep!(particles_list[1], particles_list[2], smooth_kernel, xsph_wc_lm, dr_forward_euler, step);
        catch e
            println(e);
            println("step: ", step);
            writeVTP(step, step * dr_forward_euler.dt_, vtp_io, particles_list);
            break;
        end
        if isOutputStep(step, dr_forward_euler)
            writeVTP(div(step, dr_forward_euler.output_step_), step * dr_forward_euler.dt_, vtp_io, particles_list);
        end
    end
    return nothing;
end