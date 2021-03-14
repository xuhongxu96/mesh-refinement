find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonDataModel
  IOGeometry
  REQUIRED
)

add_library(algo STATIC
    judgers/RefineJudger.h
    judgers/RefineJudger.cpp

    interpolaters/IDW.h
    interpolaters/IDW.cpp

    refiners/IRefiner.h
    refiners/InterpolateRefiner.h
    refiners/InterpolateRefiner.cpp
    refiners/DelaunayRefiner.h
    refiners/DelaunayRefiner.cpp

    optimizers/AngleOptimizer.h
    optimizers/AngleOptimizer.cpp
    optimizers/ConnectionOptimizer.h
    optimizers/ConnectionOptimizer.cpp

    io/MeshReader.h
    io/MeshReader.cpp
    io/MeshWriter.h
    io/MeshWriter.cpp
    io/GRDReader.h
    io/GRDReader.cpp
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
