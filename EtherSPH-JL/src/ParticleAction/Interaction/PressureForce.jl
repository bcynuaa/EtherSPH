#=
  @ author: bcynuaa
  @ date: 2023-11-25 15:42:13
  @ description:
 =#

function pressureForce!(
    p_i::AbstractParticle,
    p_j::AbstractParticle,
    particle_neighbour::ParticleNeighbour,
    equation_model::AbstractEquationModel
)::Nothing
    return nothing;
end

function pressureForce!(
    p_i::FluidParticle,
    p_j::FluidParticle,
    particle_neighbour::ParticleNeighbour,
    wc_liquid_model::WeakCompressibleLiquidModel
)::Nothing
    pressure_force_vec::AbstractVector = 
        -(p_i.p_ / p_i.rho_^2 + p_j.p_ / p_j.rho_^2) * particle_neighbour.kernel_gradient_vec_;
    p_i.dv_vec_ .+= p_j.mass_ * pressure_force_vec;
    p_j.dv_vec_ .-= p_i.mass_ * pressure_force_vec;
    return nothing;
end

function pressureForce!(
    p_i::FluidParticle,
    p_j::WallParticle,
    particle_neighbour::ParticleNeighbour,
    wc_liquid_model::WeakCompressibleLiquidModel
)::Nothing
    if p_i.rho_ < wc_liquid_model.rho_0_
        return nothing;
    else
        p_i.dv_vec_ .+= -2*p_i.p_ * p_i.mass_ / p_i.rho_^2 .* particle_neighbour.kernel_gradient_vec_;
        return nothing;
    end
end

function pressureForce!(
    p_i::WallParticle,
    p_j::FluidParticle,
    particle_neighbour::ParticleNeighbour,
    wc_liquid_model::WeakCompressibleLiquidModel
)::Nothing
    pressureForce!(p_j, p_i, -particle_neighbour, wc_liquid_model);
    return nothing;
end