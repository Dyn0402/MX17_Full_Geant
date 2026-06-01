cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(MX17FullSim VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)

add_executable(mx17_full_sim mx17_full_sim.cc ${SOURCES})
target_link_libraries(mx17_full_sim ${Geant4_LIBRARIES})

if(ROOT_FOUND)
    target_include_directories(mx17_full_sim PRIVATE ${ROOT_INCLUDE_DIRS})
    target_link_libraries(mx17_full_sim ${ROOT_LIBRARIES})
    target_compile_definitions(mx17_full_sim PRIVATE USE_ROOT)
    message(STATUS "ROOT found: ${ROOT_VERSION}")
else()
    message(STATUS "ROOT not found — CSV output only")
endif()

file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})

install(TARGETS mx17_full_sim DESTINATION bin)
