
cmake_minimum_required(VERSION 3.14)
project(GFoil)

set(CMAKE_CXX_STANDARD 17)
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

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.12.0
)
FetchContent_MakeAvailable(pybind11)

#############################################
# Unified physics headers (include/)
#############################################
set(UNIFIED_INCLUDES ${PROJECT_SOURCE_DIR}/include)
set(NOISE_INCLUDES   ${PROJECT_SOURCE_DIR}/noise_includes)

set(COMMON_COMPILE_DEFS
  -DUSE_CODIPACK
  -DCODI_EnableEigen
  -DDO_BL_GRADIENT=${DO_BL_GRADIENT_DEFINE}
)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(COMMON_OPTS -O2 -fno-fast-math)
endif()

#############################################
# Shared utility sources (non-template, non-main)
#############################################
set(UTIL_SOURCES
  ${PROJECT_SOURCE_DIR}/srcAD/error_package.cpp
  ${PROJECT_SOURCE_DIR}/srcAD/fresenel.cpp
)

#############################################
# Forward solver executable (GFoil_fwd_codi)
# Uses unified include/, Real = codi::RealReverse
#############################################
add_executable(GFoil_fwd_codi src/main.cpp ${UTIL_SOURCES})
target_link_libraries(GFoil_fwd_codi PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)
target_include_directories(GFoil_fwd_codi PRIVATE
  ${UNIFIED_INCLUDES}
  ${NOISE_INCLUDES}
  ${codipack_SOURCE_DIR}/include
)
target_compile_definitions(GFoil_fwd_codi PRIVATE ${COMMON_COMPILE_DEFS} -DFWD_CODI_VERSION)
target_compile_options(GFoil_fwd_codi PRIVATE ${COMMON_OPTS})

#############################################
# Adjoint solver executable (GFoil_AD)
# Uses unified include/, Real = codi::RealReverseVec<N>
#############################################
add_executable(GFoil_AD srcAD/main.cpp ${UTIL_SOURCES})
target_link_libraries(GFoil_AD PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)
target_include_directories(GFoil_AD PRIVATE
  ${UNIFIED_INCLUDES}
  ${NOISE_INCLUDES}
  ${codipack_SOURCE_DIR}/include
)
target_compile_definitions(GFoil_AD PRIVATE ${COMMON_COMPILE_DEFS} -DAD_VERSION)
target_compile_options(GFoil_AD PRIVATE ${COMMON_OPTS})

#############################################
# Python extension module (gfoil_core)
# Direct pybind11 bindings — no subprocess/JSON overhead
#############################################
pybind11_add_module(gfoil_core python/gfoil_bindings.cpp ${UTIL_SOURCES})
target_link_libraries(gfoil_core PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)
target_include_directories(gfoil_core PRIVATE
  ${UNIFIED_INCLUDES}
  ${NOISE_INCLUDES}
  ${codipack_SOURCE_DIR}/include
)
target_compile_definitions(gfoil_core PRIVATE ${COMMON_COMPILE_DEFS})
target_compile_options(gfoil_core PRIVATE ${COMMON_OPTS})
