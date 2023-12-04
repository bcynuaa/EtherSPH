#=
  @ author: bcynuaa
  @ date: 2023-11-24 16:55:12
  @ description:
 =#

function updateDensity!(p::AbstractParticle, dt::Double)::Nothing
    return nothing;
end

function updateDensity!(p::MovableParticle, dt::Double)::Nothing
    p.rho_ += p.drho_ * dt;
    p.drho_ = 0.;
    return nothing;
end