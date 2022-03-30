cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(3D_rendering)
SET(VTK_DIR ADD PATH TO VTK BUILD HERE)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(3D_rendering MACOSX_BUNDLE 3D_rendering)

if(VTK_LIBRARIES)
  target_link_libraries(3D_rendering ${VTK_LIBRARIES})
else()
  target_link_libraries(3D_rendering vtkHybrid)
endif()

