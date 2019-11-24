cmake_minimum_required(VERSION 3.7.2)
project(HMMLIB)



#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")


 
include_directories(
       ${CMAKE_SOURCE_DIR}/src/
       ${CMAKE_SOURCE_DIR}/src/Decoding
       ${CMAKE_SOURCE_DIR}/src/Evaluation
       ${CMAKE_SOURCE_DIR}/src/Evaluation/blas
       ${CMAKE_SOURCE_DIR}/src/Evaluation/conventional
       ${CMAKE_SOURCE_DIR}/src/Evaluation/csr
       ${CMAKE_SOURCE_DIR}/src/Evaluation/sblas
       ${CMAKE_SOURCE_DIR}/src/HMM
       ${CMAKE_SOURCE_DIR}/src/Training
       /usr/local/include)

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

# Initialize shared object
add_library(HMMLIB SHARED ${PROJECT_SOURCE_DIR}/src/shared.c)

target_link_libraries(HMMLIB m)
target_link_libraries(HMMLIB cblas)
 
# RSB dependency

target_link_libraries(HMMLIB rsb)
#target_link_libraries(HMMLIB hwloc)
#target_link_libraries(HMMLIB z)
target_link_libraries(HMMLIB gfortran)


target_link_libraries(HMMLIB /usr/local/lib/librsb.a)
#target_link_libraries(HMMLIB /usr/local/include)






set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE 1) #or #set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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
       ${CMAKE_SOURCE_DIR}/src/Training
       /usr/local/include)


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
target_link_libraries(exeHMMLIB m cblas)
target_link_libraries(exeHMMLIB rsb)
#target_link_libraries(exeHMMLIB hwloc)
#target_link_libraries(exeHMMLIB z)
target_link_libraries(exeHMMLIB gfortran)


target_link_libraries(exeHMMLIB /usr/local/lib/librsb.a)
#target_link_libraries(exeHMMLIB /usr/local/include)







#
#
# For codecoverage 
#
#
if(CMAKE_BUILD_TYPE MATCHES Debug)
       SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
       SET(GCC_COVERAGE_LINK_FLAGS    "-lgcov")
endif(CMAKE_BUILD_TYPE MATCHES Debug)


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





