

cmake_minimum_required(VERSION 3.14)
project(gradFoil)

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

add_executable(CFoil_fwd_codi src/main.cpp ${SOURCES})
target_link_libraries(CFoil_fwd_codi PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)

target_include_directories(CFoil_fwd_codi PRIVATE
 ${COMMON_INCLUDES_FWD}
 ${NOISE_INCLUDES_FWD}
 ${codipack_SOURCE_DIR}/include
)

target_compile_definitions(CFoil_fwd_codi PRIVATE
  -DUSE_CODIPACK
  -DCODI_EnableEigen
  -DFWD_CODI_VERSION
  -DDO_BL_GRADIENT=${DO_BL_GRADIENT_DEFINE}
)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  target_compile_options(CFoil_fwd_codi PRIVATE 
    -O2
    -fno-fast-math
  )
endif()

######## the AD version ######

set(COMMON_INCLUDES_AD ${PROJECT_SOURCE_DIR}/srcAD/include)
set(NOISE_INCLUDES_AD ${PROJECT_SOURCE_DIR}/srcAD/noise_includes)
add_executable(CFoil_AD srcAD/main.cpp srcAD/error_package.cpp)
target_link_libraries(CFoil_AD PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)
target_include_directories(CFoil_AD PRIVATE
 ${COMMON_INCLUDES_AD}
 ${NOISE_INCLUDES_AD}
 ${codipack_SOURCE_DIR}/include
)
target_compile_definitions(CFoil_AD PRIVATE
  -DUSE_CODIPACK
  -DCODI_EnableEigen
  -DAD_VERSION
  -DDO_BL_GRADIENT=${DO_BL_GRADIENT_DEFINE}
)
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  target_compile_options(CFoil_AD PRIVATE 
    -O2
    -fno-fast-math
  )
endif()