#=
  @ author: bcynuaa
  @ date: 2024-01-05 17:41:41
  @ description:
 =#

# Hagen Poiseuille flow in 2D

include("../../src/EtherSPH.jl")

const pipe_ratio = 10;
const pipe_diameter = 0.01;
const pipe_length = pipe_ratio * pipe_diameter;

const influence_ratio = 3;
const dr = pipe_diameter / 50;
const influence_radius = dr * influence_ratio;
const dim = 2;

const wall_thick_number = influence_ratio;
const buffer_thick_number = influence_ratio;

const fluid_length_number = pipe_length / dr |> round |> Int;
const fluid_width_number = pipe_diameter / dr |> round |> Int;

smooth_kernel = SmoothKernel(influence_radius, dim, CubicSplineKernel);

const ax = 8e-5;
const body_force_vec = [ax, 0.];
const rho_0 = 1000.;
const c_0 = 0.008;
const p_0 = 0.;
const mu_0 = 1e-3;
const gamma = 7.;
const wc_lm = CommonWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, body_force_vec);

const u_max = ax / 8 / mu_0 * rho_0 * pipe_diameter^2;

const dt = 0.1 * smooth_kernel.influence_radius_ / c_0;
const total_time = 17.0;
const total_step = total_time / dt |> round |> Int;
const output_step = 10;
const density_reinitialized_step = 5;
dr_forward_euler = DensityReinitializedForwardEuler(dt, total_step, output_step, density_reinitialized_step);

const step_digit = div(total_step, output_step) |> Int |> string |> length;
const file_name = "poiseuille_2d";
const file_suffix = ".vtp";
const dir_path = "./Poiseuille2DData/";
const fileld_symbols = [:rho_, :p_, :c_];
const fileld_names = ["Density", "Pressure", "SoundSpeed"];
vtp_io = VTPIO(step_digit, file_name, file_suffix, dir_path, fileld_symbols, fileld_names);

const x0 = 0.;
const y0 = 0.;
const fluid_number = fluid_length_number * fluid_width_number;

function createRectangleParticles(particle_type, dx, x_0, y_0, width_number, height_number)
    particles = [particle_type(2) for i in 1: width_number * height_number];
    for i_row in 1: height_number, j_col in 1: width_number
        index = (i_row - 1) * width_number + j_col;
        particles[index].x_vec_ = [x_0 + (j_col - 1) * dx, y_0 + (i_row - 1) * dx] .+ dx / 2;
    end
    return particles;
end

fluid_particles = createRectangleParticles(LiquidParticle, dr, x0, y0, fluid_length_number, fluid_width_number);
for index in eachindex(fluid_particles)
    fluid_particles[index].rho_ = rho_0;
    updatePressure!(fluid_particles[index], wc_lm);
    fluid_particles[index].c_ = c_0;
    fluid_particles[index].mass_ = fluid_particles[index].rho_ * dr^2;
end

# bottom wall particles
bottom_wall_particles = createRectangleParticles(
    VelocityParticle, dr,
    x0,
    y0 - wall_thick_number * dr,
    fluid_length_number,
    wall_thick_number
);
for index in eachindex(bottom_wall_particles)
    bottom_wall_particles[index].normal_vec_ = [0., 1.];
    bottom_wall_particles[index].gap_ = dr;
end

# top wall particles
top_wall_particles = createRectangleParticles(
    VelocityParticle, dr,
    x0,
    y0 + fluid_width_number * dr,
    fluid_length_number,
    wall_thick_number
);
for index in eachindex(top_wall_particles)
    top_wall_particles[index].normal_vec_ = [0., -1.];
    top_wall_particles[index].gap_ = dr;
end

wall_particles = vcat(bottom_wall_particles, top_wall_particles);

# inlet particles
inlet_particles = createRectangleParticles(
    LiquidParticle, dr,
    x0 - wall_thick_number * dr,
    y0,
    wall_thick_number,
    fluid_width_number
);
for index in eachindex(inlet_particles)
    inlet_particles[index].rho_ = rho_0;
    updatePressure!(inlet_particles[index], wc_lm);
    inlet_particles[index].c_ = c_0;
    inlet_particles[index].mass_ = inlet_particles[index].rho_ * dr^2;
end

# outlet particles
outlet_particles = createRectangleParticles(
    LiquidParticle, dr,
    x0 + fluid_length_number * dr,
    y0,
    wall_thick_number,
    fluid_width_number
);
for index in eachindex(outlet_particles)
    outlet_particles[index].rho_ = rho_0;
    updatePressure!(outlet_particles[index], wc_lm);
    outlet_particles[index].c_ = c_0;
    outlet_particles[index].mass_ = outlet_particles[index].rho_ * dr^2;
end

function check(p::LiquidParticle)::Bool
    if p.x_vec_[1] > x0 + (fluid_length_number + buffer_thick_number) * dr
        return false;
    elseif p.x_vec_[1] < x0
        return false;
    else
        return true;
    end
end

function outOfInletBuffer(p::LiquidParticle)::Bool
    return p.x_vec_[1] > x0;
end

function modifyInletBuffer(p::LiquidParticle)::Nothing
    p.x_vec_[1] -= influence_radius;
    if norm(p.v_vec_) > u_max
        p.v_vec_ = [u_max, 0.];
    end
    return nothing;
end

function intoOutletBuffer(p::LiquidParticle)::Bool
    return p.x_vec_[1] > x0 + fluid_length_number * dr;
end

function outOfOutletBuffer(p::LiquidParticle)::Bool
    return p.x_vec_[1] > x0 + (fluid_length_number + buffer_thick_number) * dr;
end

main() = @inbounds @fastmath begin
    solve!(
        fluid_particles,
        wall_particles,
        inlet_particles,
        outlet_particles,
        smooth_kernel,
        wc_lm,
        dr_forward_euler,
        vtp_io,
        check,
        outOfInletBuffer,
        modifyInletBuffer,
        intoOutletBuffer,
        outOfOutletBuffer
    );
end