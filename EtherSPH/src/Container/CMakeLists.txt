cmake_minimum_required(VERSION 3.5)
project(${CONTAINER_NAME_})

aux_source_directory(./ SRC)

include_directories(${INCLUDE_PATH_})

set(LIBRARY_OUTPUT_PATH ${LIBRARY_PATH_})

add_library(${CONTAINER_LIBRARY_NAME_} SHARED ${SRC})