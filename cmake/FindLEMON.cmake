# find_path( LEMON_INCLUDE_DIRS NAMES min_cost_arborescence.h HINTS
# "$ENV{LEMON_HOME}/include" "/home/matheusota/programs/lemon/include")

find_path(
  LEMON_INCLUDE_DIRS
  NAMES lemon/dfs.h
  HINTS "/home/mjota/lemon/include")

find_library(
  LEMON_LIBRARY
  NAMES libemon.a
  HINTS "/home/mjota/lemon/lib" "$ENV{LEMON_HOME}/lib" "/home/matheusota/programs/lemon/lib")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LEMON DEFAULT_MSG LEMON_LIBRARY)
