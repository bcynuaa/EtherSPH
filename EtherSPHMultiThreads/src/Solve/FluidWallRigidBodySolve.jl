#=
  @ author: bcynuaa
  @ date: 2024-02-08 18:28:32
  @ description:
 =#

 function eachStep!(
    fluid_particles::FluidParticlesType where FluidParticlesType <: AbstractVector{<:FluidParticle},
    wall_particles::WallParticlesType where WallParticlesType <: AbstractVector{<:WallParticle},
    rigid_body::RigidBodyType where RigidBodyType <: AbstractRigidBody,
    rigid_body_particles::RigidBodyParticlesType where RigidBodyParticlesType <: AbstractVector{<:RigidBodyParticle},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModel,
    fti_dr_forward_euler::FixedDensityReinitializedForwardEuler,
    check::Function,
    fluid_neighbour_system::InPlaceNeighborList,
    fluid_wall_neighbour_system::InPlaceNeighborList,
    fluid_rigid_body_neighbour_system::InPlaceNeighborList,
    rigid_body_wall_neighbour_system::InPlaceNeighborList,
    step::IntType where IntType <: Integer
)::Nothing
    update!(fluid_neighbour_system, [p.x_vec_ for p in fluid_particles]; cutoff=smooth_kernel.influence_radius_);
    update!(fluid_wall_neighbour_system, [p.x_vec_ for p in fluid_particles], [p.x_vec_ for p in wall_particles]; cutoff=smooth_kernel.influence_radius_);
    update!(fluid_rigid_body_neighbour_system, [p.x_vec_ for p in fluid_particles], [p.x_vec_ for p in rigid_body_particles]; cutoff=smooth_kernel.influence_radius_);
    update!(rigid_body_wall_neighbour_system, [p.x_vec_ for p in rigid_body_particles], [p.x_vec_ for p in wall_particles]; cutoff=smooth_kernel.influence_radius_);
    fluid_neighbours = findNeighbours!(
        neighborlist!(fluid_neighbour_system),
        fluid_particles, fluid_particles,
        smooth_kernel
    );
    fluid_wall_neighbours = findNeighbours!(
        neighborlist!(fluid_wall_neighbour_system),
        fluid_particles, wall_particles,
        smooth_kernel
    );
    fluid_rigid_body_neighbours = findNeighbours!(
        neighborlist!(fluid_rigid_body_neighbour_system),
        fluid_particles, rigid_body_particles,
        smooth_kernel
    );
    rigid_body_wall_neighbours = findNeighbours!(
        neighborlist!(rigid_body_wall_neighbour_system),
        rigid_body_particles, wall_particles,
        smooth_kernel
    );
    Threads.@threads for f_neighbour in fluid_neighbours
        continuity!(fluid_particles[f_neighbour.i_], fluid_particles[f_neighbour.j_], f_neighbour);
    end
    Threads.@threads for f_p in fluid_particles
        updateDensity!(f_p, fti_dr_forward_euler.dt_);
        updatePressure!(f_p, wc_lm);
    end
    Threads.@threads for f_neighbour in fluid_neighbours
        pressureForce!(fluid_particles[f_neighbour.i_], fluid_particles[f_neighbour.j_], f_neighbour, smooth_kernel);
        viscosityForce!(fluid_particles[f_neighbour.i_], fluid_particles[f_neighbour.j_], f_neighbour, smooth_kernel, wc_lm);
    end
    Threads.@threads for f_w_neighbour in fluid_wall_neighbours
        wallForce!(fluid_particles[f_w_neighbour.i_], wall_particles[f_w_neighbour.j_], f_w_neighbour, smooth_kernel);
        viscosityForce!(fluid_particles[f_w_neighbour.i_], wall_particles[f_w_neighbour.j_], f_w_neighbour, smooth_kernel, wc_lm);
    end
    Threads.@threads for f_r_b_neighbour in fluid_rigid_body_neighbours
        wallForce!(fluid_particles[f_r_b_neighbour.i_], rigid_body_particles[f_r_b_neighbour.j_], f_r_b_neighbour, smooth_kernel);
        pressureForce!(fluid_particles[f_r_b_neighbour.i_], rigid_body_particles[f_r_b_neighbour.j_], f_r_b_neighbour, smooth_kernel, wc_lm);
        viscosityForce!(fluid_particles[f_r_b_neighbour.i_], rigid_body_particles[f_r_b_neighbour.j_], f_r_b_neighbour, smooth_kernel, wc_lm);
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
    Threads.@threads for r_b_neighbour in rigid_body_wall_neighbours
        wallForce!(rigid_body_particles[r_b_neighbour.i_], wall_particles[r_b_neighbour.j_], r_b_neighbour, smooth_kernel);
    end
    joinForce!(rigid_body, rigid_body_particles);
    applyBodyForce!(rigid_body, wc_lm.body_force_vec_);
    rigidBodyMotion!(rigid_body, rigid_body_particles, fti_dr_forward_euler.dt_);
    return nothing;
end

function solve!(
    fluid_particles::FluidParticlesType where FluidParticlesType <: AbstractVector{<:FluidParticle},
    wall_particles::WallParticlesType where WallParticlesType <: AbstractVector{<:WallParticle},
    rigid_body::RigidBodyType where RigidBodyType <: AbstractRigidBody,
    rigid_body_particles::RigidBodyParticlesType where RigidBodyParticlesType <: AbstractVector{<:RigidBodyParticle},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModel,
    fti_dr_forward_euler::FixedDensityReinitializedForwardEuler,
    vtp_io::VTPIO,
    check::Function
)::Nothing
    assureDirPathExist(vtp_io);
    writeVTP(0, 0., vtp_io, [fluid_particles, wall_particles, rigid_body_particles]);
    fluid_neighbour_system = InPlaceNeighborList(
        x=[p.x_vec_ for p in fluid_particles], 
        cutoff=smooth_kernel.influence_radius_,
        parallel=true
    );
    fluid_wall_neighbour_system = InPlaceNeighborList(
        x=[p.x_vec_ for p in fluid_particles], 
        y=[p.x_vec_ for p in wall_particles],
        cutoff=smooth_kernel.influence_radius_,
        parallel=true
    );
    fluid_rigid_body_neighbour_system = InPlaceNeighborList(
        x=[p.x_vec_ for p in fluid_particles], 
        y=[p.x_vec_ for p in rigid_body_particles],
        cutoff=smooth_kernel.influence_radius_,
        parallel=true
    );
    rigid_body_wall_neighbour_system = InPlaceNeighborList(
        x=[p.x_vec_ for p in rigid_body_particles], 
        y=[p.x_vec_ for p in wall_particles],
        cutoff=smooth_kernel.influence_radius_,
        parallel=true
    );
    for step in ProgressBar(1: fti_dr_forward_euler.total_step_)
        eachStep!(
            fluid_particles, wall_particles, rigid_body, rigid_body_particles,
            smooth_kernel, wc_lm, fti_dr_forward_euler, check,
            fluid_neighbour_system, fluid_wall_neighbour_system,
            fluid_rigid_body_neighbour_system, rigid_body_wall_neighbour_system,
            step
        );
        # try
        #     eachStep!(
        #         fluid_particles, wall_particles, rigid_body, rigid_body_particles,
        #         smooth_kernel, wc_lm, fti_dr_forward_euler, check,
        #         fluid_neighbour_system, fluid_wall_neighbour_system,
        #         fluid_rigid_body_neighbour_system, rigid_body_wall_neighbour_system,
        #         step
        #     );
        # catch e
        #     println(e);
        #     println("Error occurs at step: ", step);
        #     writeVTP(step, step*fti_dr_forward_euler.dt_, vtp_io, [fluid_particles, wall_particles, rigid_body_particles]);
        #     break;
        # end
        if isOutputStep(step, fti_dr_forward_euler)
            writeVTP(
                div(step, fti_dr_forward_euler.output_step_), step*fti_dr_forward_euler.dt_, 
                vtp_io, [fluid_particles, wall_particles, rigid_body_particles]
            );
        end
    end
    return nothing;
end