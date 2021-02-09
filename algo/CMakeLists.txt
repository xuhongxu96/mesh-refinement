find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonDataModel
  IOGeometry
  REQUIRED
)

add_library(algo STATIC
    RefineJudger.h
    RefineJudger.cpp

    interpolaters/IDW.h
    interpolaters/IDW.cpp

    refiners/IRefiner.h
    refiners/InterpolateRefiner.h
    refiners/InterpolateRefiner.cpp
    refiners/DelaunayRefiner.h
    refiners/DelaunayRefiner.cpp

    io/MeshReader.h
    io/MeshReader.cpp
    io/MeshWriter.h
    io/MeshWriter.cpp
    io/GRDWriter.h
    io/GRDWriter.cpp
)

target_include_directories(algo PUBLIC .)
target_link_libraries(algo PUBLIC ${VTK_LIBRARIES})
target_compile_features(algo PUBLIC cxx_std_17)

vtk_module_autoinit(
    TARGETS algo
    MODULES ${VTK_LIBRARIES}
)
