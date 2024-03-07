#=
  @ author: bcynuaa
  @ date: 2024-02-01 16:10:49
  @ description:
 =#

function eachStep!(
    fluid_particles::FluidParticlesType where FluidParticlesType <: AbstractVector{<:FluidParticle},
    velocity_particles::VelocityParticlesType where VelocityParticlesType <: AbstractVector{<:VelocityParticle},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModel,
    fti_dr_forward_euler::FixedDensityReinitializedForwardEuler,
    check::Function,
    fluid_neighbour_system::InPlaceNeighborList,
    fluid_velocity_neighbour_system::InPlaceNeighborList,
    step::IntType where IntType <: Integer
)::Nothing
    # update neighbour list
    update!(fluid_neighbour_system, [p.x_vec_ for p in fluid_particles]; cutoff=smooth_kernel.influence_radius_);
    update!(fluid_velocity_neighbour_system, [p.x_vec_ for p in fluid_particles], [p.x_vec_ for p in velocity_particles]; cutoff=smooth_kernel.influence_radius_);
    fluid_neighbours = findNeighbours!(
        neighborlist!(fluid_neighbour_system),
        fluid_particles, fluid_particles,
        smooth_kernel
    );
    fluid_velocity_neighbours = findNeighbours!(
        neighborlist!(fluid_velocity_neighbour_system),
        fluid_particles, velocity_particles,
        smooth_kernel
    );
    @floop for f_neighbour in fluid_neighbours
        continuity!(fluid_particles[f_neighbour.i_], fluid_particles[f_neighbour.j_], f_neighbour);
    end
    @floop for f_p in fluid_particles
        updateDensity!(f_p, fti_dr_forward_euler.dt_);
        updatePressure!(f_p, wc_lm);
    end
    @floop for f_neighbour in fluid_neighbours
        pressureForce!(fluid_particles[f_neighbour.i_], fluid_particles[f_neighbour.j_], f_neighbour, smooth_kernel);
        viscosityForce!(fluid_particles[f_neighbour.i_], fluid_particles[f_neighbour.j_], f_neighbour, smooth_kernel, wc_lm);
    end
    @floop for f_v_neighbour in fluid_velocity_neighbours
        wallForce!(fluid_particles[f_v_neighbour.i_], velocity_particles[f_v_neighbour.j_], f_v_neighbour, smooth_kernel);
        pressureForce!(fluid_particles[f_v_neighbour.i_], f_v_neighbour, smooth_kernel, wc_lm);
        viscosityForce!(fluid_particles[f_v_neighbour.i_], velocity_particles[f_v_neighbour.j_], f_v_neighbour, smooth_kernel, wc_lm);
    end
    @floop for f_p in fluid_particles
        updateVelocity!(f_p, fti_dr_forward_euler.dt_, wc_lm.body_force_vec_);
        updatePosition!(f_p, fti_dr_forward_euler.dt_);
    end
    if isDensityReinitializedStep(step, fti_dr_forward_euler)
        reconstructScalar!(fluid_particles, :rho_, fluid_neighbours, smooth_kernel);
    end
    if isChekStep(step, fti_dr_forward_euler)
        out_of_bounds_index = Int64[];
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
    wc_lm::WeaklyCompressibleLiquidModel,
    fti_dr_forward_euler::FixedDensityReinitializedForwardEuler,
    vtp_io::VTPIO,
    check::Function
)::Nothing
    assureDirPathExist(vtp_io);
    writeVTP(0, 0., vtp_io, [fluid_particles, velocity_particles]);
    fluid_neighbour_system = InPlaceNeighborList(
        x=[p.x_vec_ for p in fluid_particles], 
        cutoff=smooth_kernel.influence_radius_,
        parallel=true
    );
    fluid_velocity_neighbour_system = InPlaceNeighborList(
        x=[p.x_vec_ for p in fluid_particles], 
        y=[p.x_vec_ for p in velocity_particles],
        cutoff=smooth_kernel.influence_radius_,
        parallel=true
    );
    for step in ProgressBar(1: fti_dr_forward_euler.total_step_)
        eachStep!(
            fluid_particles, velocity_particles, 
            smooth_kernel, wc_lm, fti_dr_forward_euler, 
            check,
            fluid_neighbour_system, fluid_velocity_neighbour_system,
            step
        );
        try
            eachStep!(
                fluid_particles, velocity_particles, 
                smooth_kernel, wc_lm, fti_dr_forward_euler, 
                check,
                fluid_neighbour_system, fluid_velocity_neighbour_system,
                step
            );
        catch e
            println(e);
            println("Error occurs at step: ", step);
            writeVTP(step, step*fti_dr_forward_euler.dt_, vtp_io, [fluid_particles, velocity_particles]);
            break;
        end
        if isOutputStep(step, fti_dr_forward_euler)
            writeVTP(
                div(step, fti_dr_forward_euler.output_step_), step*fti_dr_forward_euler.dt_, 
                vtp_io, [fluid_particles, velocity_particles]
            );
        end
    end
    return nothing;
end