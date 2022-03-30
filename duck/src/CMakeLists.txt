cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(duck)
SET(VTK_DIR ADD PATH TO VTK BUILD HERE)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(duck MACOSX_BUNDLE duck)

if(VTK_LIBRARIES)
  target_link_libraries(duck ${VTK_LIBRARIES})
else()
  target_link_libraries(duck vtkHybrid)
endif()

