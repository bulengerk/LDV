cmake_minimum_required(VERSION 3.15)
project(gwosc_plotter LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(CURL REQUIRED)

add_executable(gwosc_plotter main.cpp)
target_include_directories(gwosc_plotter PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(gwosc_plotter PRIVATE ${HDF5_LIBRARIES} CURL::libcurl)
if (UNIX)
  target_link_libraries(gwosc_plotter PRIVATE m)
endif()

