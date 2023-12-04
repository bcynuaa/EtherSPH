cmake_minimum_required(VERSION 3.5)
project(${PARTICLE_NAME_})

aux_source_directory(./ SRC)

include_directories(${INCLUDE_PATH_})

set(LIBRARY_OUTPUT_PATH ${LIBRARY_PATH_})

add_library(${PARTICLE_LIBRARY_NAME_} SHARED ${SRC})
target_link_libraries(${PARTICLE_NAME_} ${CONST_VARIAVLE_NAME_})
target_link_libraries(${PARTICLE_NAME_} ${CONTAINER_NAME_})
target_link_libraries(${PARTICLE_NAME_} ${KERNEL_FUNCTION_NAME_})