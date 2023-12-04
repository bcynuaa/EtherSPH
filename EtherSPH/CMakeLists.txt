# cmake minimum version
cmake_minimum_required(VERSION 3.5)

# project name
project(Ether-SPH VERSION 0.1.0 LANGUAGES CXX)

# set C++ standard as 11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# definitions
add_definitions(-DDEBUG_MODE_)

# set forlder and module name
set(CONST_VARIABLE_NAME_ ConstVariable)
set(CONTAINER_NAME_ Container)
set(KERNEL_FUNCTION_NAME_ KernelFunction)
set(PARTICLE_NAME_ Particle)
set(SOLVER_NAME_ Solver)

# set library name
set(CONST_VARIABLE_LIBRARY_NAME_ ${CONST_VARIABLE_NAME_})
set(CONTAINER_LIBRARY_NAME_ ${CONTAINER_NAME_})
set(KERNEL_FUNCTION_LIBRARY_NAME_ ${KERNEL_FUNCTION_NAME_})
set(PARTICLE_LIBRARY_NAME_ ${PARTICLE_NAME_})
set(SOLVER_LIBRARY_NAME_ ${SOLVER_NAME_})

# source path
set(SOURCE_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/src)
# include path
set(INCLUDE_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/include)
# library path
set(LIBRARY_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/lib)
# test path
set(TEST_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/test)
# draft path
set(DRAFT_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/draft)
# executable path
set(EXECUTABLE_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# add subdirectory for src to build library
add_subdirectory(${SOURCE_PATH_}/${CONST_VARIABLE_NAME_})
add_subdirectory(${SOURCE_PATH_}/${CONTAINER_NAME_})
add_subdirectory(${SOURCE_PATH_}/${KERNEL_FUNCTION_NAME_})
add_subdirectory(${SOURCE_PATH_}/${PARTICLE_NAME_})
add_subdirectory(${SOURCE_PATH_}/${SOLVER_NAME_})

# add subdirectory for test to build executable
add_subdirectory(${TEST_PATH_}/${CONTAINER_NAME_}/Array)
add_subdirectory(${TEST_PATH_}/${CONTAINER_NAME_}/NumArray)

add_subdirectory(${TEST_PATH_}/${KERNEL_FUNCTION_NAME_}/CubicSpline)
add_subdirectory(${TEST_PATH_}/${KERNEL_FUNCTION_NAME_}/Gaussian)
add_subdirectory(${TEST_PATH_}/${KERNEL_FUNCTION_NAME_}/Wendland)
add_subdirectory(${TEST_PATH_}/${KERNEL_FUNCTION_NAME_}/WendlandC4)

add_subdirectory(${TEST_PATH_}/${PARTICLE_NAME_}/ParticleStruct)
add_subdirectory(${TEST_PATH_}/${PARTICLE_NAME_}/SelfAction)
add_subdirectory(${TEST_PATH_}/${PARTICLE_NAME_}/InteractiveAction)

add_subdirectory(${TEST_PATH_}/${SOLVER_NAME_}/WaterWallSolver)