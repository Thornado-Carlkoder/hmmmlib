cmake_minimum_required(VERSION 3.7.2)
project(HMMLIB)

include_directories(
       ${CMAKE_SOURCE_DIR}/src/
       ${CMAKE_SOURCE_DIR}/src/Decoding
       ${CMAKE_SOURCE_DIR}/src/Evaluation
       ${CMAKE_SOURCE_DIR}/src/Evaluation/blas
       ${CMAKE_SOURCE_DIR}/src/Evaluation/conventional
       ${CMAKE_SOURCE_DIR}/src/Evaluation/csr
       ${CMAKE_SOURCE_DIR}/src/Evaluation/sblas
       ${CMAKE_SOURCE_DIR}/src/HMM
       ${CMAKE_SOURCE_DIR}/src/Training)

file(GLOB SOURCES
       ${PROJECT_SOURCE_DIR}/src/shared.c
       ${PROJECT_SOURCE_DIR}/src/Decoding/*.c
       ${PROJECT_SOURCE_DIR}/src/Decoding/*.h
       ${PROJECT_SOURCE_DIR}/src/Evaluation/blas/*.c
       ${PROJECT_SOURCE_DIR}/src/Evaluation/blas/*.h
       ${PROJECT_SOURCE_DIR}/src/Evaluation/conventional/*.c
       ${PROJECT_SOURCE_DIR}/src/Evaluation/conventional/*.h
       ${PROJECT_SOURCE_DIR}/src/Evaluation/csr/*.c
       ${PROJECT_SOURCE_DIR}/src/Evaluation/csr/*.h
       ${PROJECT_SOURCE_DIR}/src/Evaluation/sblas/*.c
       ${PROJECT_SOURCE_DIR}/src/Evaluation/sblas/*.h
       ${PROJECT_SOURCE_DIR}/src/HMM/*.c
       ${PROJECT_SOURCE_DIR}/src/HMM/*.h
       ${PROJECT_SOURCE_DIR}/src/Training/*.c
       ${PROJECT_SOURCE_DIR}/src/Training/*.h)

add_library(HMMLIB SHARED ${PROJECT_SOURCE_DIR}/src/shared.c)
target_link_libraries(HMMLIB m)
target_link_libraries(HMMLIB cblas)

set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE 1) #or #set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#
#
# For codecoverage
#
#
#SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")

include_directories(
       ${CMAKE_SOURCE_DIR}/src/
       ${CMAKE_SOURCE_DIR}/src/Decoding
       ${CMAKE_SOURCE_DIR}/src/Evaluation
       ${CMAKE_SOURCE_DIR}/src/Evaluation/blas
       ${CMAKE_SOURCE_DIR}/src/Evaluation/conventional
       ${CMAKE_SOURCE_DIR}/src/Evaluation/csr
       ${CMAKE_SOURCE_DIR}/src/Evaluation/sblas
       ${CMAKE_SOURCE_DIR}/src/HMM
       ${CMAKE_SOURCE_DIR}/Test
       ${CMAKE_SOURCE_DIR}/src/Training)


file(GLOB exeSOURCES
       ${PROJECT_SOURCE_DIR}/src/main.c
       ${PROJECT_SOURCE_DIR}/src/Decoding/*.c
       ${PROJECT_SOURCE_DIR}/src/Decoding/*.h
       ${PROJECT_SOURCE_DIR}/src/Evaluation/blas/*.c
       ${PROJECT_SOURCE_DIR}/src/Evaluation/blas/*.h
       ${PROJECT_SOURCE_DIR}/src/Evaluation/conventional/*.c
       ${PROJECT_SOURCE_DIR}/src/Evaluation/conventional/*.h
       ${PROJECT_SOURCE_DIR}/src/Evaluation/csr/*.c
       ${PROJECT_SOURCE_DIR}/src/Evaluation/csr/*.h
       ${PROJECT_SOURCE_DIR}/src/Evaluation/sblas/*.c
       ${PROJECT_SOURCE_DIR}/src/Evaluation/sblas/*.h
       ${PROJECT_SOURCE_DIR}/src/HMM/*.c
       ${PROJECT_SOURCE_DIR}/src/HMM/*.h
       ${PROJECT_SOURCE_DIR}/Test/*.c
       ${PROJECT_SOURCE_DIR}/Test/*.h
       ${PROJECT_SOURCE_DIR}/src/Training/*.c
       ${PROJECT_SOURCE_DIR}/src/Training/*.h)

add_executable(exeHMMLIB ${exeSOURCES})
set_source_files_properties(exeHMMLIB PROPERTIES COMPILE_FLAGS "-fprofile-arcs -ftest-coverage" )
target_link_libraries(exeHMMLIB m cblas)


#
#
# Creates a filestructure to be used by the IDE.
#
#
foreach(FILE ${exeSOURCES}) 
    # Get the directory of the source file
    get_filename_component(PARENT_DIR "${FILE}" DIRECTORY)

    # Remove common directory prefix to make the group
    string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}" "" GROUP "${PARENT_DIR}")

    # Make sure we are using windows slashes
    string(REPLACE "/" "\\" GROUP "${GROUP}")

    # Group into "Source Files" and "Header Files"
    if ("${FILE}" MATCHES ".*\\.cpp")
       set(GROUP "Source Files${GROUP}")
    elseif("${FILE}" MATCHES ".*\\.h")
       set(GROUP "Header Files${GROUP}")
    endif()

    source_group("${GROUP}" FILES "${FILE}")
endforeach()





