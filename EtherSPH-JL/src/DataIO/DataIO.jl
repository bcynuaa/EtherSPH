#=
  @ author: bcynuaa
  @ date: 2023-12-06 16:51:23
  @ description:
 =#

using WriteVTK;
import Dates;

const wall_time_format = "yyyy/mm/dd HH:MM:SS.SSSS";

abstract type AbstractDataIO end;

function assureDirPathExist(data_io::DataIOType where DataIOType<:AbstractDataIO)::Nothing
    if !isdir(data_io.dir_path_)
        mkdir(data_io.dir_path_);
    end
    return nothing;
end

include("./VTPIO.jl");