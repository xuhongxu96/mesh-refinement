find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonDataModel
  IOGeometry
  REQUIRED
)

add_library(algo STATIC
    MeshReader.h
    MeshReader.cpp
    Refiner.h
    Refiner.cpp
)

target_include_directories(algo PUBLIC .)
target_link_libraries(algo PUBLIC ${VTK_LIBRARIES})
target_compile_features(algo PUBLIC cxx_std_17)

vtk_module_autoinit(
    TARGETS algo
    MODULES ${VTK_LIBRARIES}
)
