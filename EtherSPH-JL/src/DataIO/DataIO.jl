#=
  @ author: bcynuaa
  @ date: 2023-11-28 15:01:21
  @ description:
 =#

abstract type AbstractDataIO end

function assureDirPathExist(
    data_io::AbstractDataIO
)::Nothing
    if !isdir(data_io.dir_path_)
        mkdir(data_io.dir_path_);
    end
    return nothing;
end

include("./VTPIO.jl"); 