#=
  @ author: bcynuaa
  @ date: 2023-11-24 16:51:11
  @ description:
 =#

function updatePosition!(p::AbstractParticle, dt::Double)::Nothing
    return nothing;
end

function updatePosition!(p::MovableParticle, dt::Double)::Nothing
    p.x_vec_ .+= (p.v_vec_ .- p.dv_vec_ .* dt / 2.) .* dt;
    p.dv_vec_ .= 0.;
    return nothing;
end