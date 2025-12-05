set(GENERAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB GENERAL_SOURCE_LIST CONFIGURE_DEPENDS "${GENERAL_DIR}/src/*.cpp")
set(GENERAL_INCLUDE_FOLDER "${GENERAL_DIR}/include")

# Note: "general" is a reserved keyword in CMake, so we use "general_lib"
add_library(general_lib STATIC ${GENERAL_SOURCE_LIST})
target_include_directories(
  general_lib PUBLIC ${GENERAL_INCLUDE_FOLDER} ${GUROBI_INCLUDE_DIRS} ${LEMON_INCLUDE_DIRS})
target_link_libraries(
  general_lib PUBLIC ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY}
                     ${LEMON_LIBRARY})
