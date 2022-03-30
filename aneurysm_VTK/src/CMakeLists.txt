cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(aneurysm_VTK)
SET(VTK_DIR ADD PATH TO VTK BUILD HERE)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(aneurysm_VTK MACOSX_BUNDLE aneurysm_VTK)

if(VTK_LIBRARIES)
  target_link_libraries(aneurysm_VTK ${VTK_LIBRARIES})
else()
  target_link_libraries(aneurysm_VTK vtkHybrid)
endif()

