#=
  @ author: bcynuaa
  @ date: 2023-11-24 16:38:51
  @ description:
 =#

function updateVelocity!(p::AbstractParticle, dt::Double)::Nothing
    return nothing;
end

function updateVelocity!(p::AbstractParticle, dt::Double, body_force_vec::AbstractVector)::Nothing
    return nothing;
end

function updateVelocity!(p::MovableParticle, dt::Double)::Nothing
    p.v_vec_ .+= p.dv_vec_ .* dt;
    return nothing;
end

function updateVelocity!(p::MovableParticle, dt::Double, body_force_vec::AbstractVector)::Nothing
    p.dv_vec_ .+= body_force_vec;
    updateVelocity!(p, dt);
    return nothing;
end