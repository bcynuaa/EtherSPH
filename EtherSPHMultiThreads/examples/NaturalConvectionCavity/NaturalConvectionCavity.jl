#=
  @ author: bcynuaa
  @ date: 2024-02-28 15:36:21
  @ description:
 =#

include("../../src/EtherSPHMultiThreads.jl"); 
using .EtherSPHMultiThreads;

const prandtl_number = 0.71;
const rayleigh_number = 1e5;
const reynolds_number = 50.;

const cavity_length = 1.;
const rho_0 = 1.;
const c_0 = 2.;
const umax = 0.05;
const nu_0 = umax * cavity_length / reynolds_number;
const mu_0 = rho_0 * nu_0;
const p_0 = 0.01 * rho_0 * c_0^2;
const gamma = 7;
const alpha = nu_0 / prandtl_number;

const g = 1.;
const g_vec = [0., -g];
const t_left = 1.;
const t_right = 0.;
const delta_t = t_left - t_right;
const t_0 = (t_left + t_right) / 2;
const kappa = 1.;
const cp = prandtl_number * kappa / mu_0;
const beta = rayleigh_number * nu_0 * alpha / g / delta_t / cavity_length^3;

const th_wc_lm = ThermalWeaklyCompressibleLiquidModel(
    rho_0, c_0, p_0, gamma, mu_0, g_vec, 
    t_0, kappa, cp, beta
);

const dim = 2;
const dr = 0.02;
const gap = dr;
const influence_radius = 3. * dr;

const smooth_kernel = SmoothKernel(influence_radius, dim, CubicSpline);

const dt = 0.1 * smooth_kernel.influence_radius_ / c_0;
const total_time = 1e2;
const total_step = total_time / dt |> round |> Int;
const output_step = 100;
const density_reinitialized_step = 5;
const check_step = 100;
const dr_forward_euler = FixedDensityReinitializedForwardEuler(
    dt, total_step, output_step, density_reinitialized_step, check_step
);

const step_digit = div(total_step, output_step) |> Int |> string |> length;
const file_name = "natural_convection_cavity";
const file_suffix = ".vtp";
const dir_path = "./NaturalConvectionCavityData/";
const fileld_symbols = [:rho_, :p_, :c_, :t_];
const field_names = ["Density", "Pressure", "SoundSpeed", "Temperature"];
vtp_io = VTPIO(step_digit, file_name, file_suffix, dir_path, fileld_symbols, field_names);

const x0 = 0.;
const y0 = 0.;

function createRectangleParticles(particle_type, dx, x_0, y_0, width_number, height_number)
    particles = [particle_type(Float64, 2) for i in 1: width_number * height_number];
    for i_row in 1: height_number, j_col in 1: width_number
        index = (i_row - 1) * width_number + j_col;
        particles[index].x_vec_ = [x_0 + (j_col - 1) * dx, y_0 + (i_row - 1) * dx] .+ dx / 2;
    end
    return particles;
end

fluid_row_number = cavity_length / dr |> round |> Int;
fluid_col_number = cavity_length / dr |> round |> Int;
fluid_particles = createRectangleParticles(LiquidThermalParticleMTh, dr, x0, y0, fluid_col_number, fluid_row_number);
for index in eachindex(fluid_particles)
    fluid_particles[index].rho_ = rho_0;
    fluid_particles[index].p_ = p_0;
    fluid_particles[index].c_ = c_0;
    fluid_particles[index].t_ = t_0;
    fluid_particles[index].mass_ = fluid_particles[index].rho_ * dr^2;
    fluid_particles[index].gap_ = gap;
    fluid_particles[index].kappa_ = kappa;
end

const cavity_thick_number = 3;
left_thermostat_particles = createRectangleParticles(
    ThermostaticWallParticle, dr,
    x0 - influence_radius,
    y0 - influence_radius,
    cavity_thick_number,
    fluid_row_number + 2 * cavity_thick_number,
);
for index in eachindex(left_thermostat_particles)
    left_thermostat_particles[index].t_ = t_left;
    left_thermostat_particles[index].kappa_ = kappa;
    left_thermostat_particles[index].gap_ = gap;
    left_thermostat_particles[index].normal_vec_ = [1., 0.];
end

right_thermostat_particles = createRectangleParticles(
    ThermostaticWallParticle, dr,
    x0 + cavity_length,
    y0 - influence_radius,
    cavity_thick_number,
    fluid_row_number + 2 * cavity_thick_number,
);
for index in eachindex(right_thermostat_particles)
    right_thermostat_particles[index].t_ = t_right;
    right_thermostat_particles[index].kappa_ = kappa;
    right_thermostat_particles[index].gap_ = gap;
    right_thermostat_particles[index].normal_vec_ = [-1., 0.];
end

bottom_wall_particles = createRectangleParticles(
    CompulsiveWallParticle, dr,
    x0,
    y0 - influence_radius,
    fluid_col_number,
    cavity_thick_number,
);
for index in eachindex(bottom_wall_particles)
    bottom_wall_particles[index].gap_ = gap;
    bottom_wall_particles[index].normal_vec_ = [0., 1.];
end

top_wall_particles = createRectangleParticles(
    CompulsiveWallParticle, dr,
    x0,
    y0 + cavity_length,
    fluid_col_number,
    cavity_thick_number,
);
for index in eachindex(top_wall_particles)
    top_wall_particles[index].gap_ = gap;
    top_wall_particles[index].normal_vec_ = [0., -1.];
end

thermostat_particles = vcat(left_thermostat_particles, right_thermostat_particles);
wall_particles = vcat(bottom_wall_particles, top_wall_particles);

check(particle::AbstractParticle)::Bool = true;

main() = @inbounds @fastmath solve!(
    fluid_particles, thermostat_particles, wall_particles, 
    smooth_kernel, th_wc_lm, dr_forward_euler, vtp_io, check
);