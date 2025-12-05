find_path(
  GUROBI_INCLUDE_DIRS
  NAMES gurobi_c++.h
  HINTS "/home/mjota/gurobi1203/linux64/include" "/home/matheus/Programs/gurobi1102/linux64/include" ${GUROBI_DIR}
        "$ENV{GUROBI_HOME}/include" "/opt/uw/gurobi/8.1.1/include")

find_library(
  GUROBI_LIBRARY
  NAMES gurobi gurobi120 gurobi110 gurobi90 gurobi81
  HINTS "/home/mjota/gurobi1203/linux64/lib" "/home/matheus/Programs/gurobi1102/linux64/lib" ${GUROBI_DIR}
        "$ENV{GUROBI_HOME}/lib" "/opt/uw/gurobi/8.1.1/lib")
find_library(
  GUROBI_CXX_LIBRARY
  NAMES gurobi_c++
  HINTS "/home/mjota/gurobi1203/linux64/lib" "/home/matheus/Programs/gurobi1102/linux64/lib" ${GUROBI_DIR}
        "$ENV{GUROBI_HOME}/lib" "/opt/uw/gurobi/8.1.1/lib")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY)
