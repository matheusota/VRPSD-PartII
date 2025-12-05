set(NODE_SCENARIOS_MODEL_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB NODE_SCENARIOS_MODEL_SOURCE_LIST CONFIGURE_DEPENDS
     "${NODE_SCENARIOS_MODEL_DIR}/src/*.cpp")
set(NODE_SCENARIOS_MODEL_INCLUDE_FOLDER "${NODE_SCENARIOS_MODEL_DIR}/include")

add_library(node_scenarios_model STATIC ${NODE_SCENARIOS_MODEL_SOURCE_LIST})
target_include_directories(
    node_scenarios_model PUBLIC ${NODE_SCENARIOS_MODEL_INCLUDE_FOLDER} ${LEMON_INCLUDE_DIRS})
target_link_libraries(
    node_scenarios_model
  PUBLIC ${NODE_SCENARIOS_LIBRARY} ${GUROBI_CXX_LIBRARY} cvrpsep
         general_lib cuts)
