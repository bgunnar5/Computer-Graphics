cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(colors)
SET(VTK_DIR ADD PATH TO VTK BUILD HERE)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(colors MACOSX_BUNDLE colors)

if(VTK_LIBRARIES)
  target_link_libraries(colors ${VTK_LIBRARIES})
else()
  target_link_libraries(colors vtkHybrid)
endif()

