#=
  @ author: bcynuaa
  @ date: 2023-11-28 20:42:05
  @ description:
 =#

include("../../src/EtherSPH.jl");

const dt = 1e-5;
const total_step = 10000;
const output_step = 100;
const dim = 2;
const dr = 1.5e-2;
const h = 3 * dr;

const g = 9.8;
const body_force_vec = [0.0, -g];
const c_0 = 50.;
const rho_0 = 1000.;
const p_0 = 0.;
const mu_0 = 1e-3;

const kernel = WendlandC2Kernel(h, dim);
const wc_liquid_model = WeakCompressibleLiquidModel(c_0, rho_0, p_0, mu_0, dim);
const forward_euler = ForwardEuler(dt, total_step, output_step);
const vtp_io = VTPIO(4, "WaterColumn", ".vtp", "./result/");

const row_fluid_particle_number = 40;
const col_fluid_particle_number = 30;
const fluid_particle_number = row_fluid_particle_number * col_fluid_particle_number;

const water_max_height = row_fluid_particle_number * 2 * dr;

fluid_particles = [LiquidParticle(dim) for i in 1: fluid_particle_number];

for i_row = 1: row_fluid_particle_number
    for j_col = 1: col_fluid_particle_number
        index = (i_row - 1) * col_fluid_particle_number + j_col;
        fluid_particles[index].x_vec_ .= [2*j_col-1, 2*i_row-1] .* dr;
        rho_g_h = rho_0 * g * (water_max_height - fluid_particles[index].x_vec_[2])
        fluid_particles[index].p_ = rho_g_h + p_0;
        fluid_particles[index].rho_ = rho_g_h / c_0^2 + rho_0;
        fluid_particles[index].mass_ = rho_0 * dr^2 * 4;
    end
end

function createBoxBoundary()
    each_edge_len_number = 80;
    each_edge_thick_number = 6;
    each_edge_number = each_edge_len_number * each_edge_thick_number;
    bottom_edge_particles = [WallParticle(dim) for i in 1: each_edge_number];
    for i in 1: each_edge_len_number, j in 1: each_edge_thick_number
        index = (i - 1) * each_edge_thick_number + j;
        bottom_edge_particles[index].x_vec_ .= [2*i-1, -2*j+1] .* dr;
        bottom_edge_particles[index].x_vec_[1] -= each_edge_thick_number * 2 * dr;
    end
    top_edge_particles = [WallParticle(dim) for i in 1: each_edge_number];
    for i in 1: each_edge_len_number, j in 1: each_edge_thick_number
        index = (i - 1) * each_edge_thick_number + j;
        top_edge_particles[index].x_vec_ .= [2*i-1, 2*j-1] .* dr;
        top_edge_particles[index].x_vec_[2] += 2*dr*(each_edge_len_number-each_edge_thick_number);
    end
    left_edge_particles = [WallParticle(dim) for i in 1: each_edge_number];
    for i in 1: each_edge_len_number, j in 1: each_edge_thick_number
        index = (i - 1) * each_edge_thick_number + j;
        left_edge_particles[index].x_vec_ .= [-2*j+1, 2*i-1] .* dr;
    end
    right_edge_particles = [WallParticle(dim) for i in 1: each_edge_number];
    for i in 1: each_edge_len_number, j in 1: each_edge_thick_number
        index = (i - 1) * each_edge_thick_number + j;
        right_edge_particles[index].x_vec_ .= [2*j-1, 2*i-1] .* dr;
        right_edge_particles[index].x_vec_[1] += 2*dr*(each_edge_len_number-each_edge_thick_number);
        right_edge_particles[index].x_vec_[2] -= 2*each_edge_thick_number*dr;
    end
    return vcat(bottom_edge_particles, top_edge_particles, left_edge_particles, right_edge_particles);
end

particle_pool = vcat(fluid_particles, createBoxBoundary()) |> ParticlePool;

main() = @inbounds @fastmath solve!(particle_pool, body_force_vec, kernel, wc_liquid_model, forward_euler, vtp_io);