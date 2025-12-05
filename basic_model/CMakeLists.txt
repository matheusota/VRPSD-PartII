set(BASIC_MODEL_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB BASIC_MODEL_SOURCE_LIST CONFIGURE_DEPENDS
     "${BASIC_MODEL_DIR}/src/*.cpp")
set(BASIC_MODEL_INCLUDE_FOLDER "${BASIC_MODEL_DIR}/include")

add_library(basic_model STATIC ${BASIC_MODEL_SOURCE_LIST})
target_include_directories(basic_model PUBLIC ${BASIC_MODEL_INCLUDE_FOLDER})
target_link_libraries(
  basic_model PUBLIC ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} cvrpsep
                     general_lib cuts)
