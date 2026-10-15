cmake_minimum_required(VERSION 3.21)
project(make_it_stand_cherchi_core LANGUAGES C CXX)

# Standalone optional backend build. No download and no experimental source
# inclusion: both the Boolean revision and the interpreter are explicit.
get_filename_component(MAKE_IT_STAND_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(CHERCHI_SOURCE "${MAKE_IT_STAND_ROOT}/external/reference-sources-20260917/sources/InteractiveAndRobustMeshBooleans-7bd6c2697695f555f9a15d7e154e084d841316f0")
if(NOT EXISTS "${CHERCHI_SOURCE}/code/booleans.h")
  message(FATAL_ERROR "Pinned Cherchi source missing: ${CHERCHI_SOURCE}")
endif()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_VARIABLE pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY)
find_package(pybind11 3.0.1 EXACT CONFIG REQUIRED)
set(TBB_TEST OFF CACHE BOOL "" FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory("${CHERCHI_SOURCE}/arrangements/external/oneTBB" "${CMAKE_BINARY_DIR}/bundled-tbb")
set(cinolib_DIR "${CHERCHI_SOURCE}/arrangements/external/Cinolib")
set(CINOLIB_USES_OPENGL_GLFW_IMGUI OFF CACHE BOOL "" FORCE)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON CACHE BOOL "" FORCE)
find_package(cinolib CONFIG REQUIRED)
pybind11_add_module(_cherchi_core_boolean MODULE NO_EXTRAS InProcessBoolean.cpp)
target_include_directories(_cherchi_core_boolean PRIVATE
  "${CHERCHI_SOURCE}" "${CHERCHI_SOURCE}/code" "${CHERCHI_SOURCE}/arrangements/code"
  "${CHERCHI_SOURCE}/arrangements/external/Indirect_Predicates/include"
  "${CHERCHI_SOURCE}/arrangements/external/abseil-cpp")
target_compile_definitions(_cherchi_core_boolean PRIVATE TBB_PARALLEL=1)
target_compile_options(_cherchi_core_boolean PRIVATE -O2 -frounding-math -mavx2)
target_link_libraries(_cherchi_core_boolean PRIVATE cinolib tbb)
set_target_properties(_cherchi_core_boolean PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python")

option(MAKE_IT_STAND_CHERCHI_RATIONAL_ARITHMETIC "Build exact GMP surface moments and transport coefficients" OFF)
if(MAKE_IT_STAND_CHERCHI_RATIONAL_ARITHMETIC)
  find_path(MAKE_IT_STAND_GMPXX_INCLUDE_DIR gmpxx.h REQUIRED)
  find_library(MAKE_IT_STAND_GMPXX_LIBRARY gmpxx REQUIRED)
  find_library(MAKE_IT_STAND_GMP_LIBRARY gmp REQUIRED)
  pybind11_add_module(_cherchi_rational_moments MODULE NO_EXTRAS RationalMoments.cpp)
  target_include_directories(_cherchi_rational_moments PRIVATE "${MAKE_IT_STAND_GMPXX_INCLUDE_DIR}")
  target_compile_options(_cherchi_rational_moments PRIVATE -O3)
  target_link_libraries(_cherchi_rational_moments PRIVATE "${MAKE_IT_STAND_GMPXX_LIBRARY}" "${MAKE_IT_STAND_GMP_LIBRARY}")
  set_target_properties(_cherchi_rational_moments PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python")
endif()
