#=
  @ author: bcynuaa
  @ date: 2024-02-06 14:20:26
  @ description:
 =#

mutable struct RigidBodyParticleMTh{
    RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType}, 
    AtomArrayType <: AbstractVector{Base.Threads.Atomic{RealType}},
} <: RigidBodyParticle
    x_vec_::ArrayType
    v_vec_::ArrayType
    dv_vec_::AtomArrayType
    mass_::RealType
    normal_vec_::ArrayType
    gap_::RealType
end

function RigidBodyParticleMTh(RealType::DataType, dim::IntType)::RigidBodyParticleMTh where IntType <: Integer
    x_vec = zeros(RealType, dim);
    v_vec = zeros(RealType, dim);
    dv_vec = [Base.Threads.Atomic{RealType}(0.) for i in 1: dim];
    mass::RealType = 0.;
    normal_vec = zeros(RealType, dim);
    gap::RealType = 0.;
    RigidBodyParticleMTh(
        x_vec, v_vec, dv_vec, 
        mass, normal_vec, gap
    );
end

abstract type AbstractRigidBody end;

mutable struct RigidBody2D{
    RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType},
    AtomArrayType <: AbstractVector{Base.Threads.Atomic{RealType}},
} <: AbstractRigidBody
    x_vec_::ArrayType # centroid position
    v_vec_::ArrayType # centroid velocity
    dv_vec_::AtomArrayType # centroid acceleration
    omega_::RealType # angular velocity
    alpha_::Base.Threads.Atomic{RealType} # angular acceleration
    mass_::RealType # mass
    inertia_::RealType # inertia
end

function RigidBody2D(RealType::DataType)::RigidBody2D
    x_vec = zeros(RealType, 2);
    v_vec = zeros(RealType, 2);
    dv_vec = [Base.Threads.Atomic{RealType}(0.) for i in 1: 2];
    omega::RealType = 0.;
    alpha = Base.Threads.Atomic{RealType}(0.);
    mass::RealType = 0.;
    inertia::RealType = 0.;
    RigidBody2D(
        x_vec, 
        v_vec, dv_vec, 
        omega, alpha,
        mass, inertia
    );
end

function joinForce!(
    rigid_body_2d::RigidBody2D, 
    rigid_body_particles::RigidBodyParticleArrayType
)::Nothing where RigidBodyParticleArrayType <: AbstractVector{<:RigidBodyParticle}
    # ! after calculation
    # ! dv_vec_ is ∑ⱼ mⱼ aⱼ
    # ! omega_ is ∑ⱼ mⱼ (rⱼ × aⱼ)
    force_vec = similar(rigid_body_2d.x_vec_);
    Threads.@threads for particle in rigid_body_particles
        for i_dim in eachindex(force_vec)
            force_vec[i_dim] = particle.dv_vec_[i_dim][] * particle.mass_;
            Threads.atomic_add!(rigid_body_2d.dv_vec_[i_dim], force_vec[i_dim]);
            particle.dv_vec_[i_dim][] = 0.;
        end
        r_vec = particle.x_vec_ .- rigid_body_2d.x_vec_;
        Threads.atomic_add!(rigid_body_2d.alpha_, r_vec[1]*force_vec[2] - r_vec[2]*force_vec[1]);
    end
    return nothing;
end

function applyBodyForce!(
    rigid_body_2d::RigidBody2D,
    body_force_vec::ArrayType where ArrayType <: AbstractVector{<:RealType} # here body_force_vec means m/s^2
)::Nothing where RealType <: AbstractFloat
    for i_dim in eachindex(rigid_body_2d.dv_vec_)
        Base.Threads.atomic_add!(rigid_body_2d.dv_vec_[i_dim], body_force_vec[i_dim] * rigid_body_2d.mass_);
    end
    return nothing;
end

function rigidBodyMotion!(
    rigid_body_2d::RigidBody2D,
    rigid_body_particles::RigidBodyParticleArrayType where RigidBodyParticleArrayType <: AbstractVector{<:RigidBodyParticle},
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    for i_dim in eachindex(rigid_body_2d.dv_vec_)
        rigid_body_2d.v_vec_[i_dim] += rigid_body_2d.dv_vec_[i_dim][] / rigid_body_2d.mass_ * dt;
        rigid_body_2d.dv_vec_[i_dim][] = 0.;
    end
    rigid_body_2d.omega_ += rigid_body_2d.alpha_[] / rigid_body_2d.inertia_ * dt;
    rigid_body_2d.alpha_[] = 0.;
    Threads.@threads for particle in rigid_body_particles
        r_vec = particle.x_vec_ .- rigid_body_2d.x_vec_;
        tau_vec = [-r_vec[2], r_vec[1]];
        particle.v_vec_ .= rigid_body_2d.v_vec_ .+ rigid_body_2d.omega_ * tau_vec;
        particle.x_vec_ .+= particle.v_vec_ * dt;
        d_theta = rigid_body_2d.omega_ * dt;
        particle.normal_vec_ = [cos(d_theta) -sin(d_theta); sin(d_theta) cos(d_theta)] * particle.normal_vec_;
    end
    rigid_body_2d.x_vec_ .+= rigid_body_2d.v_vec_ * dt;
    return nothing;
end