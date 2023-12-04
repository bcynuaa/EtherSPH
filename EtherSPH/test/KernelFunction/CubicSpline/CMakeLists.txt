cmake_minimum_required(VERSION 3.5)
project(CubicSpline_test)

aux_source_directory(. SRC)
set(EXECUTABLE_OUTPUT_PATH ${EXECUTABLE_PATH_}/test/${KERNEL_FUNCTION_NAME_})

include_directories(${INCLUDE_PATH_})
link_directories(${LIBRARY_PATH_})

add_executable(${PROJECT_NAME} ${SRC})
target_link_libraries(${PROJECT_NAME} ${CONST_VARIAVLE_NAME_})
target_link_libraries(${PROJECT_NAME} ${CONTAINER_LIBRARY_NAME_})
target_link_libraries(${PROJECT_NAME} ${KERNEL_FUNCTION_LIBRARY_NAME_})