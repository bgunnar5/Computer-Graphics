cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(triangles)
SET(VTK_DIR ADD PATH TO VTK BUILD HERE)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(triangles MACOSX_BUNDLE triangles)

if(VTK_LIBRARIES)
  target_link_libraries(triangles ${VTK_LIBRARIES})
else()
  target_link_libraries(triangles vtkHybrid)
endif()

