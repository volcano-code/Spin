cmake_minimum_required(VERSION 3.21)
project(cherchi_reference_probe LANGUAGES C CXX)

set(CHERCHI_SOURCE "${CMAKE_CURRENT_LIST_DIR}/../../external/reference-sources-20260917/sources/InteractiveAndRobustMeshBooleans-7bd6c2697695f555f9a15d7e154e084d841316f0")
get_filename_component(CHERCHI_SOURCE "${CHERCHI_SOURCE}" ABSOLUTE)
if(NOT EXISTS "${CHERCHI_SOURCE}/code/booleans.h")
  message(FATAL_ERROR "Pinned Cherchi source missing: ${CHERCHI_SOURCE}")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(TBB_TEST OFF CACHE BOOL "" FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory("${CHERCHI_SOURCE}/arrangements/external/oneTBB" "${CMAKE_BINARY_DIR}/bundled-tbb")

set(cinolib_DIR "${CHERCHI_SOURCE}/arrangements/external/Cinolib")
set(CINOLIB_USES_OPENGL_GLFW_IMGUI OFF CACHE BOOL "" FORCE)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON CACHE BOOL "" FORCE)
find_package(cinolib CONFIG REQUIRED)

add_executable(cherchi_probe probe.cpp)
target_include_directories(cherchi_probe PRIVATE
  "${CHERCHI_SOURCE}"
  "${CHERCHI_SOURCE}/code"
  "${CHERCHI_SOURCE}/arrangements/code"
  "${CHERCHI_SOURCE}/arrangements/external/Indirect_Predicates/include"
  "${CHERCHI_SOURCE}/arrangements/external/abseil-cpp"
)
target_compile_definitions(cherchi_probe PRIVATE TBB_PARALLEL=1)
target_compile_options(cherchi_probe PRIVATE -O2 -frounding-math -mavx2)
target_link_libraries(cherchi_probe PRIVATE cinolib tbb)
