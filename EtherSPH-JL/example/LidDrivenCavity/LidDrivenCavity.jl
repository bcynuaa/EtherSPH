#=
  @ author: bcynuaa
  @ date: 2024-01-03 15:26:16
  @ description:
 =#

# Ghia et al 1980

include("../../src/EtherSPH.jl");

const reynolds_number = 100.;

const lid_length = 1.;
const lid_velocity = 1.;
const rho_0 = 1.;
const mu_0 = lid_velocity * lid_length * rho_0 / reynolds_number;
const c_0 = 10. * lid_velocity;
const p_0 = 0.;
const gamma = 7.;

const dr = 0.02;
const gap = dr;
const influence_radius = 3. * dr;
const dim = 2;

smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2Kernel);

const body_force_vec = [0., 0.];

const wc_lm = CommonWeaklyCompressibleLiquidModel(rho_0, c_0, p_0, gamma, mu_0, body_force_vec);

const dt = 0.1 * smooth_kernel.influence_radius_ / c_0;
const total_time = 15.;
const total_step = total_time / dt |> round |> Int;
const output_step = 100;
const density_reinitialized_step = 5;

dr_forward_euler = DensityReinitializedForwardEuler(dt, total_step, output_step, density_reinitialized_step);

const step_digit = div(total_step, output_step) |> Int |> string |> length;
const file_name = "lid_driven_cavity";
const file_suffix = ".vtp";
const dir_path = "./LidDrivenCavityData/";
const fileld_symbols = [:rho_, :p_, :c_];
const fileld_names = ["Density", "Pressure", "SoundSpeed"];
vtp_io = VTPIO(step_digit, file_name, file_suffix, dir_path, fileld_symbols, fileld_names);

const x0 = 0.;
const y0 = 0.;

function createRectangleParticles(particle_type, dx, x_0, y_0, width_number, height_number)
    particles = [particle_type(2) for i in 1: width_number * height_number];
    for i_row in 1: height_number, j_col in 1: width_number
        index = (i_row - 1) * width_number + j_col;
        particles[index].x_vec_ = [x_0 + (j_col - 1) * dx, y_0 + (i_row - 1) * dx] .+ dx / 2;
    end
    return particles;
end

fluid_row_number = lid_length / dr |> round |> Int;
fluid_col_number = lid_length / dr |> round |> Int;
fluid_particles = createRectangleParticles(LiquidParticle, dr, x0, y0, fluid_col_number, fluid_row_number);
for index in eachindex(fluid_particles)
    fluid_particles[index].rho_ = rho_0;
    fluid_particles[index].p_ = p_0;
    fluid_particles[index].c_ = c_0;
    fluid_particles[index].mass_ = fluid_particles[index].rho_ * dr^2;
    fluid_particles[index].gap_ = gap;
end

const lid_thick_number = 3;
const lid_length_number = lid_length / dr |> round |> Int;

lid_particles = createRectangleParticles(VelocityParticle, dr, x0, y0 + lid_length, lid_length_number, lid_thick_number);
for index in eachindex(lid_particles)
    lid_particles[index].v_vec_ = [lid_velocity, 0.];
    lid_particles[index].gap_ = gap;
    lid_particles[index].normal_vec_ = [0., -1.];
end

# left VelocityParticle
left_particles = createRectangleParticles(VelocityParticle, dr, x0 - lid_thick_number * dr, y0, lid_thick_number, lid_length_number);
for index in eachindex(left_particles)
    left_particles[index].v_vec_ = [0., 0.];
    left_particles[index].gap_ = gap;
    left_particles[index].normal_vec_ = [1., 0.];
end

# right VelocityParticle
right_particles = createRectangleParticles(VelocityParticle, dr, x0 + lid_length, y0, lid_thick_number, lid_length_number);
for index in eachindex(right_particles)
    right_particles[index].v_vec_ = [0., 0.];
    right_particles[index].gap_ = gap;
    right_particles[index].normal_vec_ = [-1., 0.];
end

# bottom VelocityParticle
bottom_particles = createRectangleParticles(VelocityParticle, dr, x0, y0 - lid_thick_number * dr, lid_length_number, lid_thick_number);
for index in eachindex(bottom_particles)
    bottom_particles[index].v_vec_ = [0., 0.];
    bottom_particles[index].gap_ = gap;
    bottom_particles[index].normal_vec_ = [0., 1.];
end

# left bottom
left_bottom_particles = createRectangleParticles(VelocityParticle, dr, x0 - lid_thick_number * dr, y0 - lid_thick_number * dr, lid_thick_number, lid_thick_number);
for index in eachindex(left_bottom_particles)
    left_bottom_particles[index].v_vec_ = [0., 0.];
    left_bottom_particles[index].gap_ = gap;
    left_bottom_particles[index].normal_vec_ = [1., 1.] ./ sqrt(2.);
end

# right bottom
right_bottom_particles = createRectangleParticles(VelocityParticle, dr, x0 + lid_length, y0 - lid_thick_number * dr, lid_thick_number, lid_thick_number);
for index in eachindex(right_bottom_particles)
    right_bottom_particles[index].v_vec_ = [0., 0.];
    right_bottom_particles[index].gap_ = gap;
    right_bottom_particles[index].normal_vec_ = [-1., 1.] ./ sqrt(2.);
end

# left top
left_top_particles = createRectangleParticles(VelocityParticle, dr, x0 - lid_thick_number * dr, y0 + lid_length, lid_thick_number, lid_thick_number);
for index in eachindex(left_top_particles)
    left_top_particles[index].v_vec_ = [lid_velocity, 0.];
    left_top_particles[index].gap_ = gap;
    left_top_particles[index].normal_vec_ = [1., -1.] ./ sqrt(2.);
end

# right top
right_top_particles = createRectangleParticles(VelocityParticle, dr, x0 + lid_length, y0 + lid_length, lid_thick_number, lid_thick_number);
for index in eachindex(right_top_particles)
    right_top_particles[index].v_vec_ = [lid_velocity, 0.];
    right_top_particles[index].gap_ = gap;
    right_top_particles[index].normal_vec_ = [-1., -1.] ./ sqrt(2.);
end

check(particle::AbstractParticle)::Bool = true;

velocity_particles = vcat(
    lid_particles, left_particles, right_particles, bottom_particles,
    left_bottom_particles, right_bottom_particles,
    left_top_particles, right_top_particles
);

main() = @inbounds @fastmath solve!(fluid_particles, velocity_particles, smooth_kernel, wc_lm, dr_forward_euler, vtp_io, check);