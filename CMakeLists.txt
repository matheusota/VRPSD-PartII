cmake_minimum_required(VERSION 3.10)

# Project name.
project(SVRP LANGUAGES CXX)

# Compiler requirements.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set build type to Debug if not specified.
message(STATUS "BUILD TYPE: ${CMAKE_BUILD_TYPE}")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

# Compiler flags.
add_compile_options(-m64 -no-pie -w -mtune=native ${ABI_FLAG})

# Enable LEMON assertions in Debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(LEMON_ENABLE_ASSERTS)
endif()

# Warnings add_compile_options(-Wall -Wextra -Wpedantic)
# add_link_options(-fsanitize=address -fno-omit-frame-pointer)

# Linker flags.
set(MY_LINK_FLAGS "-lm -lpthread -lboost_system ")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MY_LINK_FLAGS}")

# Debug and Release specific flags.
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

# External packages libraries.
list(APPEND CMAKE_MODULE_PATH "${SVRP_SOURCE_DIR}/cmake")
find_package(GUROBI REQUIRED)
find_package(LEMON REQUIRED)

set(Boost_USE_STATIC_LIBS ON)
find_package(
  Boost
  COMPONENTS regex system timer chrono
  REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

message(STATUS "GUROBI_INCLUDE_DIRS: ${GUROBI_INCLUDE_DIRS}")
message(STATUS "GUROBI_LIBRARY: ${GUROBI_LIBRARY}")
message(STATUS "GUROBI_CXX_LIBRARY: ${GUROBI_CXX_LIBRARY}")

# Add subdirectories.
add_subdirectory(CVRPSEP)
add_subdirectory(general)
add_subdirectory(cuts)
add_subdirectory(node_recourse_model)
add_subdirectory(node_scenarios_model)

# Main source and header.
set(MAIN_HEADER "main.h")
set(MAIN_SOURCE "main.cpp")

include_directories(${GUROBI_INCLUDE_DIRS})
add_executable(main ${MAIN_SOURCE} ${MAIN_HEADER})
target_link_libraries(
  main
  PRIVATE general_lib
          cuts
          node_recourse_model
          node_scenarios_model
          ${LEMON_LIBRARY}
          ${GUROBI_LIBRARY}
          ${GUROBI_CXX_LIBRARY}
          ${Boost_LIBRARIES})
