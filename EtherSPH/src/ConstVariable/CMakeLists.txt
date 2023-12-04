cmake_minimum_required(VERSION 3.5)
project(${CONST_VARIABLE_NAME_})

aux_source_directory(./ SRC)

include_directories(${INCLUDE_PATH_})

set(LIBRARY_OUTPUT_PATH ${LIBRARY_PATH_})

add_library(${CONST_VARIABLE_LIBRARY_NAME_} SHARED ${SRC})