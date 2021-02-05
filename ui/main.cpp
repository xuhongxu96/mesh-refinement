#include <MeshReader.h>
#include <vtkInteractorStyleTerrain.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

int main() {
  mr::MeshReader reader;
  mr::Mesh mesh = reader.Read(DATA_PATH "original.mesh");

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputData(mesh.poly_data);

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->SetDiffuse(0.6);
  actor->GetProperty()->SetDiffuseColor(80, 0, 0);
  actor->GetProperty()->SetSpecular(0.8);
  actor->GetProperty()->SetSpecularPower(20.0);

  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("ReadSTL");

  vtkNew<vtkInteractorStyleTerrain> style;
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->SetInteractorStyle(style);

  renderer->AddActor(actor);
  renderer->SetBackground(255, 255, 255);

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}