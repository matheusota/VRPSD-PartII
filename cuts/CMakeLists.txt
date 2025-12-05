set(CUTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB CUTS_SOURCE_LIST CONFIGURE_DEPENDS "${CUTS_DIR}/src/*.cpp")
set(CUTS_INCLUDE_FOLDER "${CUTS_DIR}/include")

add_library(cuts STATIC ${CUTS_SOURCE_LIST})
target_include_directories(
  cuts PUBLIC ${CUTS_INCLUDE_FOLDER} ${CVRPSEP_INCLUDE_FOLDER} ${LEMON_INCLUDE_DIRS})
target_link_libraries(cuts PUBLIC ${LEMON_LIBRARY} ${GUROBI_LIBRARY}
                                  ${GUROBI_CXX_LIBRARY} cvrpsep general_lib)
