#=
  @ author: bcynuaa
  @ date: 2023-12-07 21:41:00
  @ description:
 =#

 include("../../src/EtherSPH.jl");

 const fluid_row_number = 200;
 const fluid_col_number = 100;

 const dim = 2;
 const dr = 1e-2;
 const gap = dr;
 const influence_radius = 1.5 * dr;

 const fluid_height = fluid_row_number * dr;
 
 smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2Kernel);
 
 const g = 9.8;
 const g_vec = [0., -g];
 const rho_0 = 1000.;
 const c_0 = 120.;
#  const c_0 = sqrt(g * dr * fluid_row_number) * 10.;
 const gamma = 7.;
 const mu_0 = 1e-3;
 xsph_wc_lm = XSPHWeaklyCompressibleLiquidModel(rho_0, c_0, gamma, mu_0, g_vec, fluid_height);
 
 const dt = 0.1 * smooth_kernel.influence_radius_ / c_0;
#  const dt = 5e-5;
 const total_time = 3.;
 const total_step = div(total_time, dt) |> Int;
 const output_step = 100;
 const density_reinitialized_step = 10;
 dr_forward_euler = DensityReinitializedForwardEuler(dt, total_step, output_step, density_reinitialized_step);
 
 const step_digit = div(total_step, output_step) |> Int |> string |> length;
 const file_name = "collapse_dry";
 const file_suffix = ".vtp";
 const dir_path = "./CollapseDryData/";
 const fileld_symbols = [:rho_, :p_, :c_];
 const fileld_names = ["Density", "Pressure", "SoundSpeed"];
 vtp_io = VTPIO(step_digit, file_name, file_suffix, dir_path, fileld_symbols, fileld_names);
 
 const x0 = 0.;
 const y0 = 0.;
 const fluid_number = fluid_row_number * fluid_col_number;
 function createFluidParticles()
     particles = [LiquidParticle(2) for i in 1: fluid_number];
     for i_row in 1: fluid_row_number, j_col in 1: fluid_col_number
         index = (i_row - 1) * fluid_col_number + j_col;
         particles[index].x_vec_ = [x0 + (j_col - 1) * dr, y0 + (i_row - 1) * dr] .+ dr / 2;
        #  particles[index].p_ = rho_0 * g * (fluid_height - particles[index].x_vec_[2]);
        #  particles[index].rho_ = (particles[index].p_ / xsph_wc_lm.b_ + 1.)^(1. / gamma) * rho_0;
         particles[index].rho_ = rho_0;
         updatePressure!(particles[index], xsph_wc_lm);
         particles[index].mass_ = rho_0 * dr^2;
     end
     return particles;
 end
 fluid_particles = createFluidParticles();
 
 const bottom_wall_row_number = 3;
 const bottom_wall_col_number = 400;
 const bottom_wall_number = bottom_wall_row_number * bottom_wall_col_number;
 function createBottomWallParticles()
     particles = [WallParticle(2) for i in 1: bottom_wall_number];
     for i_row in 1: bottom_wall_row_number, j_col in 1: bottom_wall_col_number
         index = (i_row - 1) * bottom_wall_col_number + j_col;
         particles[index].x_vec_ = [x0 + (j_col - 1) * dr, y0 - dr * i_row] .+ dr / 2;
         particles[index].normal_vec_ = [0., 1.];
         particles[index].gap_ = dr;
     end
     return particles;
 end
 bottom_wall_particles = createBottomWallParticles();
 
 const left_wall_row_number = 300 + bottom_wall_row_number * 1;
 const left_wall_col_number = bottom_wall_row_number;
 const left_wall_number = left_wall_row_number * left_wall_col_number;
 function createLeftWallParticles()
     particles = [WallParticle(2) for i in 1: left_wall_number];
     for i_row in 1: left_wall_row_number, j_col in 1: left_wall_col_number
         index = (i_row - 1) * left_wall_col_number + j_col;
         particles[index].x_vec_ = [x0 - dr * j_col, y0 + (i_row - 1) * dr] .+ dr / 2;
         particles[index].x_vec_[2] -= dr * bottom_wall_row_number;
         particles[index].normal_vec_ = [1., 0.];
         particles[index].gap_ = dr;
     end
     for particle in particles
         if particle.x_vec_[2] < y0
            #  particle.normal_vec_ = [1., 1.] ./ sqrt(2.);
            particle.normal_vec_ = [0., 1.];
        #  elseif particle.x_vec_[2] > y0 + (left_wall_row_number - 2*bottom_wall_row_number) * dr
        #      particle.normal_vec_ = [1., -1.] ./ sqrt(2.);
         end
     end
     return particles;
 end
 left_wall_particles = createLeftWallParticles();
 
 right_wall_particles = deepcopy(left_wall_particles);
 for particle in right_wall_particles
     particle.x_vec_[1] += dr * (bottom_wall_col_number + bottom_wall_row_number);
     particle.normal_vec_ = [-1., 0.];
     if particle.x_vec_[2] < y0
        #  particle.normal_vec_ = [-1., 1.] ./ sqrt(2.);
        particle.normal_vec_ = [0., 1.];
    #  elseif particle.x_vec_[2] > y0 + (left_wall_row_number - 2*bottom_wall_row_number) * dr
    #      particle.normal_vec_ = [-1., -1.] ./ sqrt(2.);
     end
 end
 
 top_wall_particles = deepcopy(bottom_wall_particles);
 for particle in top_wall_particles
     particle.x_vec_[2] += dr * (left_wall_row_number - bottom_wall_row_number);
     particle.normal_vec_ = [0., -1.];
 end
 
 wall_particles = vcat(bottom_wall_particles, left_wall_particles, right_wall_particles);
 
 particles = [fluid_particles, wall_particles];
 
 main() = @inbounds @fastmath solve!(particles, smooth_kernel, xsph_wc_lm, dr_forward_euler, vtp_io);