project (ui)

find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonDataModel
  IOGeometry
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
  REQUIRED
)

add_executable(ui 
	main.cpp
)

target_link_libraries(ui PRIVATE algo ${VTK_LIBRARIES})
target_compile_definitions(ui PRIVATE "DATA_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/../data/\"")
vtk_module_autoinit(
    TARGETS ui
    MODULES ${VTK_LIBRARIES}
)
