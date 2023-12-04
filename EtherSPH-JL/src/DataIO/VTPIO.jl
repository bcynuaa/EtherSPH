#=
  @ author: bcynuaa
  @ date: 2023-11-28 15:01:47
  @ description:
 =#

# * use LightXML to write particles to vtp file

struct VTPIO <: AbstractDataIO
    step_digit_::Isize
    file_name_::String
    file_suffix_::String
    dir_path_::String
end

function vtpFileNameAtStep(step::Isize, vtp_io::VTPIO)::String
    return joinpath(
        vtp_io.dir_path_,
        string(
            vtp_io.file_name_,
            string(step, pad=vtp_io.step_digit_),
            vtp_io.file_suffix_
        )
    );
end

function writeVtp(
    particle_pool::ParticlePool,
    step::Isize,
    vtp_io::VTPIO
)::Nothing
    dim::Isize = length(particle_pool.particles_[1].x_vec_);
    float_type::DataType = typeof(particle_pool.particles_[1].x_vec_[1]);
    # get data from particle pool by different particle type
    points::AbstractArray{float_type} = zeros(float_type, dim, length(particle_pool.particles_));
    velocitys::AbstractArray{float_type} = zeros(float_type, dim, length(particle_pool.particles_));
    rhos::AbstractArray{float_type} = zeros(float_type, length(particle_pool.particles_));
    pressures::AbstractArray{float_type} = zeros(float_type, length(particle_pool.particles_));
    current_type_first_index::Isize = 1;
    for particle_type in keys(particle_pool.type_index_dict_)
        n_particles = length(particle_pool.type_index_dict_[particle_type]);
        particle_index = particle_pool.type_index_dict_[particle_type]
        if particle_type <: FixedParticle
            for i = 1: n_particles
                points[:, current_type_first_index + i - 1] = particle_pool.particles_[particle_index[i]].x_vec_;
            end
            # other property set as NaN
            velocitys[:, current_type_first_index:current_type_first_index + n_particles - 1] .= NaN;
            rhos[current_type_first_index:current_type_first_index + n_particles - 1] .= NaN;
            pressures[current_type_first_index:current_type_first_index + n_particles - 1] .= NaN;
        elseif particle_type <: MovableParticle
            for i = 1: n_particles
                points[:, current_type_first_index + i - 1] = particle_pool.particles_[particle_index[i]].x_vec_;
                velocitys[:, current_type_first_index + i - 1] = particle_pool.particles_[particle_index[i]].v_vec_;
                rhos[current_type_first_index + i - 1] = particle_pool.particles_[particle_index[i]].rho_;
                pressures[current_type_first_index + i - 1] = particle_pool.particles_[particle_index[i]].p_;
            end
        end
        current_type_first_index += n_particles;
    end
    # write to vtp file
    vtp_file_name::String = vtpFileNameAtStep(step, vtp_io);
    cells = [MeshCell(PolyData.Verts(), [i]) for i in 1: length(particle_pool.particles_)];
    vtk_file = vtk_grid(vtp_file_name, points, cells);
    vtk_file["Velocity"] = velocitys;
    vtk_file["Density"] = rhos;
    vtk_file["Pressure"] = pressures;
    vtk_save(vtk_file);
    return nothing;
end