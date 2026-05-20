

cmake_minimum_required(VERSION 3.14)
project(GFoil)

set(CMAKE_CXX_STANDARD 17)

# Default build type
set(CMAKE_BUILD_TYPE Release)
include(FetchContent)

option(DO_BL_GRADIENT "Enable boundary layer gradient calculation" OFF)

if(DO_BL_GRADIENT)
  set(DO_BL_GRADIENT_DEFINE 1)
else()
  set(DO_BL_GRADIENT_DEFINE 0)
endif()
#############################################
# Dependencies
#############################################
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)

FetchContent_Declare(
  codipack
  GIT_REPOSITORY https://github.com/SciCompKL/CoDiPack.git
  GIT_TAG master
)
FetchContent_MakeAvailable(codipack)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz)
FetchContent_MakeAvailable(json)

#############################################

set(COMMON_INCLUDES_FWD ${PROJECT_SOURCE_DIR}/src/include)
set(NOISE_INCLUDES_FWD ${PROJECT_SOURCE_DIR}/src/noise_includes)

file(GLOB_RECURSE SOURCES src/*.cpp)
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/main.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/gfoil_fwd_bindings.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/gfoil_ad_bindings.cpp")

add_executable(GFoil_fwd_codi src/main.cpp ${SOURCES})
target_link_libraries(GFoil_fwd_codi PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)

target_include_directories(GFoil_fwd_codi PRIVATE
 ${COMMON_INCLUDES_FWD}
 ${NOISE_INCLUDES_FWD}
 ${codipack_SOURCE_DIR}/include
)

target_compile_definitions(GFoil_fwd_codi PRIVATE
  -DUSE_CODIPACK
  -DCODI_EnableEigen
  -DFWD_CODI_VERSION
  -DDO_BL_GRADIENT=${DO_BL_GRADIENT_DEFINE}
)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  target_compile_options(GFoil_fwd_codi PRIVATE 
    -O2
    -fno-fast-math
  )
endif()

######## the AD version ######

set(COMMON_INCLUDES_AD ${PROJECT_SOURCE_DIR}/srcAD/include)
set(NOISE_INCLUDES_AD ${PROJECT_SOURCE_DIR}/srcAD/noise_includes)
add_executable(GFoil_AD srcAD/main.cpp src/error_package.cpp)
target_link_libraries(GFoil_AD PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)
target_include_directories(GFoil_AD PRIVATE
 ${COMMON_INCLUDES_AD}
 ${COMMON_INCLUDES_FWD}
 ${NOISE_INCLUDES_AD}
 ${NOISE_INCLUDES_FWD}
 ${codipack_SOURCE_DIR}/include
)
target_compile_definitions(GFoil_AD PRIVATE
  -DUSE_CODIPACK
  -DCODI_EnableEigen
  -DAD_VERSION
  -DDO_BL_GRADIENT=${DO_BL_GRADIENT_DEFINE}
)
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  target_compile_options(GFoil_AD PRIVATE
    -O2
    -fno-fast-math
  )
endif()

######## pybind11 module ########

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# All src/ sources except main.cpp (run_forward.cpp, coupled.cpp, etc.)
# SOURCES is already defined above.

pybind11_add_module(gfoil_cpp
    src/gfoil_fwd_bindings.cpp
    src/gfoil_ad_bindings.cpp
    ${SOURCES}
    src/error_package.cpp
)

target_link_libraries(gfoil_cpp PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)

target_include_directories(gfoil_cpp PRIVATE
    ${COMMON_INCLUDES_FWD}
    ${NOISE_INCLUDES_FWD}
    ${PROJECT_SOURCE_DIR}/srcAD/include
    ${PROJECT_SOURCE_DIR}/srcAD/noise_includes
    ${codipack_SOURCE_DIR}/include
)

target_compile_definitions(gfoil_cpp PRIVATE
    -DUSE_CODIPACK
    -DCODI_EnableEigen
    -DFWD_CODI_VERSION
    -DDO_BL_GRADIENT=${DO_BL_GRADIENT_DEFINE}
)

target_compile_options(gfoil_cpp PRIVATE -O2 -fno-fast-math)

# Place the .so in the Python package directory so `from . import gfoil_cpp` works
set_target_properties(gfoil_cpp PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/GFoil"
)