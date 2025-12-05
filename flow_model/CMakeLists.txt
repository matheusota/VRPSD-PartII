set(FLOW_MODEL_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB FLOW_MODEL_SOURCE_LIST CONFIGURE_DEPENDS
     "${FLOW_MODEL_DIR}/src/*.cpp")
set(FLOW_MODEL_INCLUDE_FOLDER "${FLOW_MODEL_DIR}/include")

add_library(flow_model STATIC ${FLOW_MODEL_SOURCE_LIST})
target_include_directories(flow_model PUBLIC ${FLOW_MODEL_INCLUDE_FOLDER})
target_link_libraries(
  flow_model PUBLIC ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} ${LEMON_LIBRARY}
                    cvrpsep general_lib cuts)
