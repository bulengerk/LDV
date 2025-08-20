cmake_minimum_required(VERSION 3.16)
project(gwosc_plotter LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Wconversion)
endif()

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(CURL REQUIRED)

add_library(gw STATIC
  src/util.cpp
  src/http.cpp
  src/h5.cpp
  src/dsp.cpp
  src/plot.cpp
  src/csv.cpp
  src/app.cpp
)

target_include_directories(gw PUBLIC include ${HDF5_INCLUDE_DIRS})
target_link_libraries(gw PUBLIC ${HDF5_LIBRARIES} CURL::libcurl)
if (UNIX)
  target_link_libraries(gw PUBLIC m)
endif()

add_executable(gwosc_plotter src/main.cpp)
target_link_libraries(gwosc_plotter PRIVATE gw)
