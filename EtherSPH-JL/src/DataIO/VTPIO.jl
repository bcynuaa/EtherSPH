#=
  @ author: bcynuaa
  @ date: 2023-12-06 16:54:33
  @ description:
 =#

struct VTPIO <: AbstractDataIO
    step_digit_::IntType where IntType<:Integer
    file_name_::String
    file_suffix_::String
    dir_path_::String
    field_symbols_::FieldSymbolsType where FieldSymbolsType <: AbstractVector{Symbol}
    field_names_::FieldNamesType where FieldNamesType <: AbstractVector{String}
end

function vtpFileNameAtStep(
    step::IntType where IntType<:Integer,
    vtp_io::VTPIO
)::String
    return joinpath(
        vtp_io.dir_path_,
        string(
            vtp_io.file_name_,
            string(step, pad=vtp_io.step_digit_),
            vtp_io.file_suffix_
        )
    );
end

function getParticleFields(
    field_symbols::FieldSymbolsType where FieldSymbolsType <: AbstractVector{Symbol},
    particle::ParticleArrayType where ParticleArrayType <: AbstractVector{<:AbstractParticle}
)
    field_number = length(field_symbols);
    particle_number = length(particle);
    fields = zeros(field_number, particle_number);
    particle_field_symbols = fieldnames(eltype(particle));
    for i_field in eachindex(field_symbols)
        if field_symbols[i_field] in particle_field_symbols
            fields[i_field, :] = [getfield(particle[i_particle], field_symbols[i_field]) for i_particle in eachindex(particle)];
        else
            fields[i_field, :] .= NaN;
        end
    end
    return fields;
end

function writeVTP(
    step::IntType,
    t::RealType,
    vtp_io::VTPIO,
    particles_list::ParticlesListType where ParticlesListType <: AbstractVector{<:AbstractVector}
)::Nothing where {IntType <: Integer, RealType <: AbstractFloat}
    dim::IntType = length(particles_list[1][1].x_vec_);
    n_particles_list = [length(particles) for particles in particles_list];
    n_particles::IntType = sum(n_particles_list);
    points = zeros(RealType, dim, n_particles);
    velocitys = zeros(RealType, dim, n_particles);
    fields = zeros(RealType, length(vtp_io.field_symbols_), n_particles);
    current_point_index::IntType = 1;
    for particles in particles_list
        particle_number = length(particles);
        fields[:, current_point_index: current_point_index + particle_number - 1] = getParticleFields(vtp_io.field_symbols_, particles);
        if eltype(particles) <: MovableParticle
            for i_particle in eachindex(particles)
                points[:, current_point_index + i_particle - 1] = particles[i_particle].x_vec_;
                velocitys[:, current_point_index + i_particle - 1] = particles[i_particle].v_vec_;
            end
        elseif eltype(particles) <: VelocityParticle
            for i_particle in eachindex(particles)
                points[:, current_point_index + i_particle - 1] = particles[i_particle].x_vec_;
                velocitys[:, current_point_index + i_particle - 1] = particles[i_particle].v_vec_;
            end
        else
            for i_particle in eachindex(particles)
                points[:, current_point_index + i_particle - 1] = particles[i_particle].x_vec_;
                velocitys[:, current_point_index + i_particle - 1] .= NaN;
            end
        end
        current_point_index += particle_number;
    end
    vtp_file_name::String = vtpFileNameAtStep(step, vtp_io);
    cells = [MeshCell(PolyData.Verts(), [i]) for i in 1: n_particles];
    vtk_file = vtk_grid(vtp_file_name, points, cells);
    vtk_file["Step"] = step;
    vtk_file["Time"] = t;
    vtk_file["WallTime"] = Dates.format(Dates.now(), wall_time_format);
    vtk_file["Velocity"] = velocitys;
    for i_field in eachindex(vtp_io.field_symbols_)
        vtk_file[vtp_io.field_names_[i_field]] = fields[i_field, :];
    end
    vtk_save(vtk_file);
    return nothing;
end