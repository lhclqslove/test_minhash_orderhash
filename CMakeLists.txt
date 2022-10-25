cmake_minimum_required(VERSION 3.12)
project(test_minhash_orderhash)

set(CMAKE_CXX_STANDARD 14)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CHECKS=${CHECKS} (Use -DCHECKS=True to enable all checks")
message(STATUS "LOG=${LOG} (Use -LOG=True to output log files for analysis")
message(STATUS "Directory: ${CMAKE_CURRENT_SOURCE_DIR}")
###############################
# boost
set( Boost_NO_SYSTEM_PATHS ON ) # do not use system boost
add_subdirectory(boost-cmake)
#openmp
find_package(OpenMP)
#xxhash
find_package(xxHash 0.7 CONFIG REQUIRED)

add_library(testLib
        STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/dnaToBits.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ReadData.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/dnaToBits.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/ReadData.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Types.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/BBHashMap.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/ReadFilter.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/BBHashMap.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Minhash_ReadFilter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/xxhash.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/OrderHashReadFilter.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/OrderHashReadFilter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/BooPHF.h
        include/Test_minhash.h include/Test_Orderminhash.h)
target_include_directories(testLib
        PUBLIC include src
        )
target_link_libraries(testLib PUBLIC OpenMP::OpenMP_CXX -DLIBBSC_OPENMP_SUPPORT)
target_link_libraries(testLib PRIVATE xxHash::xxhash)
target_link_libraries(testLib
        PUBLIC
            Boost::filesystem
            Boost::iostreams
            Boost::program_options
            -lz -ldl -lm -lpthread
)


#############
add_executable(test_minhash_orderhash src/main.cpp )
target_link_libraries(test_minhash_orderhash PUBLIC xxHash::xxhash)
target_link_libraries(test_minhash_orderhash PUBLIC testLib)