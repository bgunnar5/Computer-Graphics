cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(camera_manipulation)
SET(VTK_DIR ADD PATH TO VTK BUILD HERE)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(camera_manipulation MACOSX_BUNDLE camera_manipulation)

if(VTK_LIBRARIES)
  target_link_libraries(camera_manipulation ${VTK_LIBRARIES})
else()
  target_link_libraries(camera_manipulation vtkHybrid)
endif()

