#=
  @ author: bcynuaa
  @ date: 2024-01-04 14:40:42
  @ description:
 =#

# M.A. Cruchaga

# 3D example

include("../../src/EtherSPH.jl");

const box_a = 0.42;
const box_b = 0.44;
const box_c = 0.228;

const water_a = 0.114;
const water_b = 0.114;
const water_c = 0.228;

const dim = 3;

const dr = water_a / 20; # 20 * 20 * 40 = 16000
const gap = dr;
const influence_radius = 3. * dr;

smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2Kernel);

const water_a_number = water_a / dr |> round |> Int;
const water_b_number = water_b / dr |> round |> Int;
const water_c_number = water_c / dr |> round |> Int;

const box_a_number = box_a / dr |> round |> Int;
const box_b_number = box_b / dr |> round |> Int;
const box_c_number = box_c / dr |> round |> Int;

const box_thick_number = 3;

const gravity = 9.81;
const g_vec = [0., -gravity, 0.];
const rho_0 = 1e3;
const c_0 = 10.;
const p_0 = 0.;
const gamma = 7.;
const mu_0 = 1e-3;
wc_lm = CommonWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, g_vec);

const dt = 0.1 * smooth_kernel.influence_radius_ / c_0;
const total_time = 1.;
const total_step = total_time / dt |> round |> Int;
const output_step = 20;
const density_reinitialized_step = 5;
dr_forward_euler = DensityReinitializedForwardEuler(dt, total_step, output_step, density_reinitialized_step);

const step_digit = div(total_step, output_step) |> Int |> string |> length;
const file_name = "cruchaga_3d";
const file_suffix = ".vtp";
const dir_path = "./Cruchaga3DData/";
const fileld_symbols = [:rho_, :p_, :c_];
const fileld_names = ["Density", "Pressure", "SoundSpeed"];
vtp_io = VTPIO(step_digit, file_name, file_suffix, dir_path, fileld_symbols, fileld_names);

const x0 = 0.;
const y0 = 0.;
const z0 = 0.;

const water_number = water_a_number * water_b_number * water_c_number;

function createBoxParticles(particle_type, dx, x_0, y_0, z_0, x_number, y_number, z_number)
    particles = [particle_type(3) for i in 1: x_number * y_number * z_number];
    for i_z in 1: z_number, i_y in 1: y_number, i_x in 1: x_number
        index = (i_z - 1) * x_number * y_number + (i_y - 1) * x_number + i_x;
        particles[index].x_vec_ = [x_0 + (i_x - 1) * dx, y_0 + (i_y - 1) * dx, z_0 + (i_z - 1) * dx] .+ dx / 2;
    end
    return particles;
end

water_particles = createBoxParticles(LiquidParticle, dr, x0, y0, z0, water_a_number, water_b_number, water_c_number);
for index in eachindex(water_particles)
    water_particles[index].p_ = rho_0 * gravity * (water_b - water_particles[index].x_vec_[2]);
    water_particles[index].rho_ = (water_particles[index].p_ / wc_lm.b_ + 1.)^(1. / gamma) * rho_0;
    # fluid_particles[index].rho_ = rho_0;
    updatePressure!(water_particles[index], wc_lm);
    water_particles[index].mass_ = water_particles[index].rho_ * dr^3;
    water_particles[index].gap_ = dr;
end

bottom_wall_particles = createBoxParticles(
    WallParticle, dr,
    x0 - box_thick_number * dr,
    y0 - box_thick_number * dr,
    z0 - box_thick_number * dr,
    box_a_number + 2 * box_thick_number,
    box_thick_number,
    box_c_number + 2 * box_thick_number
);
for index in eachindex(bottom_wall_particles)
    bottom_wall_particles[index].normal_vec_ = [0., 1., 0.];
    bottom_wall_particles[index].gap_ = dr;
end

top_wall_particles = createBoxParticles(
    WallParticle, dr,
    x0 - box_thick_number * dr,
    y0 + box_b_number * dr,
    z0 - box_thick_number * dr,
    box_a_number + 2 * box_thick_number,
    box_thick_number,
    box_c_number + 2 * box_thick_number
);
for index in eachindex(top_wall_particles)
    top_wall_particles[index].normal_vec_ = [0., -1., 0.];
    top_wall_particles[index].gap_ = dr;
end

front_wall_particles = createBoxParticles(
    WallParticle, dr,
    x0 - box_thick_number * dr,
    y0,
    z0 - box_thick_number * dr,
    box_a_number + 2 * box_thick_number,
    box_b_number,
    box_thick_number
);
for index in eachindex(front_wall_particles)
    front_wall_particles[index].normal_vec_ = [0., 0., 1.];
    front_wall_particles[index].gap_ = dr;
end

back_wall_particles = createBoxParticles(
    WallParticle, dr,
    x0 - box_thick_number * dr,
    y0,
    z0 + box_c_number * dr,
    box_a_number + 2 * box_thick_number,
    box_b_number,
    box_thick_number
);
for index in eachindex(back_wall_particles)
    back_wall_particles[index].normal_vec_ = [0., 0., -1.];
    back_wall_particles[index].gap_ = dr;
end

left_wall_particles = createBoxParticles(
    WallParticle, dr,
    x0 - box_thick_number * dr,
    y0,
    z0,
    box_thick_number,
    box_b_number,
    box_c_number
);
for index in eachindex(left_wall_particles)
    left_wall_particles[index].normal_vec_ = [1., 0., 0.];
    left_wall_particles[index].gap_ = dr;
end

right_wall_particles = createBoxParticles(
    WallParticle, dr,
    x0 + box_a_number * dr,
    y0,
    z0,
    box_thick_number,
    box_b_number,
    box_c_number
);
for index in eachindex(right_wall_particles)
    right_wall_particles[index].normal_vec_ = [-1., 0., 0.];
    right_wall_particles[index].gap_ = dr;
end

wall_particles = vcat(
    bottom_wall_particles,
    top_wall_particles,
    front_wall_particles,
    back_wall_particles,
    left_wall_particles,
    right_wall_particles
);

check(particle::AbstractParticle)::Bool = true;

main() = @inbounds @fastmath solve!(water_particles, wall_particles, smooth_kernel, wc_lm, dr_forward_euler, vtp_io, check);