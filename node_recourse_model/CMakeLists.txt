set(NODE_RECOURSE_MODEL_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB NODE_RECOURSE_MODEL_SOURCE_LIST CONFIGURE_DEPENDS
     "${NODE_RECOURSE_MODEL_DIR}/src/*.cpp")
set(NODE_RECOURSE_MODEL_INCLUDE_FOLDER "${NODE_RECOURSE_MODEL_DIR}/include")

add_library(node_recourse_model STATIC ${NODE_RECOURSE_MODEL_SOURCE_LIST})
target_include_directories(
  node_recourse_model PUBLIC ${NODE_RECOURSE_MODEL_INCLUDE_FOLDER} ${LEMON_INCLUDE_DIRS})
target_link_libraries(
  node_recourse_model
  PUBLIC ${NODE_RECOURSE_LIBRARY} ${GUROBI_CXX_LIBRARY} cvrpsep
         general_lib cuts)
