#=
  @ author: bcynuaa
  @ date: 2024-01-05 20:12:04
  @ description:
 =#

#=
  @ author: bcynuaa
  @ date: 2024-01-03 16:06:26
  @ description:
 =#

function eachStep!(
    fluid_particles::FluidParticlesType where FluidParticlesType <: AbstractVector{<:FluidParticle},
    wall_particles::WallParticlesType where WallParticlesType <: AbstractVector{<:FixedParticle},
    inlet_particles::InletParticlesType where InletParticlesType <: AbstractVector{<:FluidParticle},
    outlet_particles::OutletParticlesType where OutletParticlesType <: AbstractVector{<:FluidParticle},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel,
    dr_forward_euler::DensityReinitializedForwardEuler,
    step::IntType where IntType <: Integer,
    check::Function,
    outOfInletBuffer::Function,
    modifyInletBuffer::Function,
    intoOutletBuffer::Function,
    outOfOutletBuffer::Function
)::Nothing
    all_f_particles = vcat(fluid_particles, inlet_particles, outlet_particles);
    all_f_neighbours = findNeighbours(all_f_particles, smooth_kernel);
    all_f_w_neighbours = findNeighbours(all_f_particles, wall_particles, smooth_kernel);
    # * continuity
    Threads.@threads for neighbour in all_f_neighbours
        continuity!(all_f_particles[neighbour.i_], all_f_particles[neighbour.j_], neighbour);
    end
    # * update density and pressure
    Threads.@threads for i_f_particle in eachindex(fluid_particles)
        updateDensity!(fluid_particles[i_f_particle], dr_forward_euler.dt_);
        updatePressure!(fluid_particles[i_f_particle], wc_lm);
    end
    # * momentum
    Threads.@threads for neighbour in all_f_neighbours
        pressureForce!(all_f_particles[neighbour.i_], all_f_particles[neighbour.j_], neighbour, smooth_kernel);
        viscosityForce!(all_f_particles[neighbour.i_], all_f_particles[neighbour.j_], neighbour, smooth_kernel, wc_lm);
    end
    Threads.@threads for neighbour in all_f_w_neighbours
        wallForce!(all_f_particles[neighbour.i_], wall_particles[neighbour.j_], neighbour, smooth_kernel);
        # viscosityForce!(all_f_particles[neighbour.i_], neighbour, smooth_kernel, wc_lm);
        viscosityForce!(all_f_particles[neighbour.i_], wall_particles[neighbour.j_], neighbour, smooth_kernel, wc_lm);
        pressureForce!(all_f_particles[neighbour.i_], neighbour, smooth_kernel, wc_lm);
    end
    # * update fluid particles
    Threads.@threads for i_f_particle in eachindex(fluid_particles)
        updateVelocity!(fluid_particles[i_f_particle], dr_forward_euler.dt_, wc_lm.body_force_vec_);
        updatePosition!(fluid_particles[i_f_particle], dr_forward_euler.dt_);
    end
    # * reinitialize density
    if isDensityReinitializedStep(step, dr_forward_euler)
        reconstructScalar!(all_f_particles, :rho_, all_f_neighbours, smooth_kernel);
        out_of_bounds_index = [];
        for i_f_particle in eachindex(fluid_particles)
            if !check(fluid_particles[i_f_particle])
                push!(out_of_bounds_index, i_f_particle);
            end
        end
        deleteat!(fluid_particles, out_of_bounds_index);
    end
    # * update inlet particles
    Threads.@threads for i_in_particle in eachindex(inlet_particles)
        inlet_particles[i_in_particle].drho_ = 0.;
        inlet_particles[i_in_particle].dv_vec_ .= 0.;
        inlet_particles[i_in_particle].rho_ = wc_lm.rho_0_;
        updateVelocity!(inlet_particles[i_in_particle], dr_forward_euler.dt_, wc_lm.body_force_vec_);
        updatePosition!(inlet_particles[i_in_particle], dr_forward_euler.dt_);
    end
    # * update outlet particles
    Threads.@threads for i_out_particle in eachindex(outlet_particles)
        # outlet_particles[i_out_particle].drho_ = 0.;
        # outlet_particles[i_out_particle].dv_vec_ .= 0.;
        updateVelocity!(outlet_particles[i_out_particle], dr_forward_euler.dt_, wc_lm.body_force_vec_);
        updatePosition!(outlet_particles[i_out_particle], dr_forward_euler.dt_);
    end
    # * modify inlet particles
    out_of_inlet_index = [];
    for i_in_particle in eachindex(inlet_particles)
        if outOfInletBuffer(inlet_particles[i_in_particle])
            push!(out_of_inlet_index, i_in_particle);
        end
    end
    for i_out_of_inlet_particle in out_of_inlet_index
        push!(fluid_particles, deepcopy(inlet_particles[i_out_of_inlet_particle]));
        modifyInletBuffer(inlet_particles[i_out_of_inlet_particle]);
    end
    # * modify inlet particles
    into_outlet_index = [];
    for i_f_particle in eachindex(fluid_particles)
        if intoOutletBuffer(fluid_particles[i_f_particle])
            push!(into_outlet_index, i_f_particle);
        end
    end
    for i_into_outlet_particle in into_outlet_index
        push!(outlet_particles, deepcopy(fluid_particles[i_into_outlet_particle]));
        outOfOutletBuffer(fluid_particles[i_into_outlet_particle]);
    end
    deleteat!(fluid_particles, into_outlet_index);
    out_of_outlet_index = [];
    for i_out_particle in eachindex(outlet_particles)
        if outOfOutletBuffer(outlet_particles[i_out_particle])
            push!(out_of_outlet_index, i_out_particle);
        end
    end
    deleteat!(outlet_particles, out_of_outlet_index);
    return nothing;
end

function solve!(
    fluid_particles::FluidParticlesType where FluidParticlesType <: AbstractVector{<:FluidParticle},
    wall_particles::WallParticlesType where WallParticlesType <: AbstractVector{<:FixedParticle},
    inlet_particles::InletParticlesType where InletParticlesType <: AbstractVector{<:FluidParticle},
    outlet_particles::OutletParticlesType where OutletParticlesType <: AbstractVector{<:FluidParticle},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel,
    dr_forward_euler::DensityReinitializedForwardEuler,
    vtp_io::VTPIO,
    check::Function,
    outOfInletBuffer::Function,
    modifyInletBuffer::Function,
    intoOutletBuffer::Function,
    outOfOutletBuffer::Function
)::Nothing
    assureDirPathExist(vtp_io);
    writeVTP(0, 0., vtp_io, [fluid_particles, wall_particles]);
    for step in ProgressBar(1: dr_forward_euler.total_step_)
        eachStep!(
            fluid_particles,
            wall_particles,
            inlet_particles,
            outlet_particles,
            smooth_kernel,
            wc_lm,
            dr_forward_euler,
            step,
            check,
            outOfInletBuffer,
            modifyInletBuffer,
            intoOutletBuffer,
            outOfOutletBuffer
        );
        # try
        #     eachStep!(
        #         fluid_particles,
        #         wall_particles,
        #         inlet_particles,
        #         outlet_particles,
        #         smooth_kernel,
        #         wc_lm,
        #         dr_forward_euler,
        #         step,
        #         check,
        #         outOfInletBuffer,
        #         modifyInletBuffer,
        #         intoOutletBuffer,
        #         outOfOutletBuffer
        #     );
        # catch e
        #     println(e);
        #     println("step: ", step);
        #     writeVTP(step, step * dr_forward_euler.dt_, vtp_io, [fluid_particles, wall_particles, inlet_particles, outlet_particles]);
        #     break;
        # end
        if isOutputStep(step, dr_forward_euler)
            writeVTP(div(step, dr_forward_euler.output_step_), step * dr_forward_euler.dt_, vtp_io, [fluid_particles, wall_particles]);
        end
    end
    return nothing;
end