#include <RefineJudger.h>
#include <interpolaters/IDW.h>
#include <io/GRDWriter.h>
#include <io/MeshReader.h>
#include <io/MeshWriter.h>
#include <refiners/DelaunayRefiner.h>
#include <refiners/InterpolateRefiner.h>
#include <vtkActor.h>
#include <vtkCellData.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSimplePointsReader.h>
#include <vtkUnsignedCharArray.h>

#include "MouseInteractorStyle.h"

void AddColor(vtkPolyData* res) {
  vtkNew<vtkUnsignedCharArray> colors;
  colors->SetNumberOfComponents(3);
  colors->SetNumberOfTuples(res->GetNumberOfCells());

  for (size_t i = 0; i < res->GetNumberOfCells(); ++i) {
    float rgb[3] = {100, 0, 0};
    colors->InsertTuple(i, rgb);
  }

  res->GetCellData()->SetScalars(colors);
}

void Render(vtkPolyData* left_mesh, vtkPolyData* right_mesh) {
  double left_viewport[4] = {0.0, 0.0, 0.5, 1.0};
  double right_viewport[4] = {0.5, 0.0, 1.0, 1.0};

  vtkNew<vtkPolyDataMapper> left_mapper;
  left_mapper->SetInputData(left_mesh);

  vtkNew<vtkPolyDataMapper> right_mapper;
  right_mapper->SetInputData(right_mesh);

  vtkNew<vtkActor> left_actor;
  left_actor->SetMapper(left_mapper);
  left_actor->GetProperty()->SetSpecular(.4);
  left_actor->GetProperty()->SetSpecularPower(30.0);

  vtkNew<vtkActor> right_actor;
  right_actor->SetMapper(right_mapper);
  right_actor->GetProperty()->SetSpecular(.4);
  right_actor->GetProperty()->SetSpecularPower(30.0);

  vtkNew<vtkRenderer> left_renderer;
  left_renderer->SetViewport(left_viewport);

  vtkNew<vtkRenderer> right_renderer;
  right_renderer->SetViewport(right_viewport);

  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(left_renderer);
  renderWindow->AddRenderer(right_renderer);
  renderWindow->SetWindowName("Mesh Refinement");

  vtkNew<MouseInteractorStyle> style;
  style->SetCurrentRenderer(right_renderer);
  style->DataMap = {
      {left_renderer, left_mesh},
      {right_renderer, right_mesh},
  };

  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->SetInteractorStyle(style);

  left_renderer->AddActor(left_actor);
  left_renderer->SetBackground(255, 255, 255);

  right_renderer->AddActor(right_actor);
  right_renderer->SetBackground(255, 255, 255);

  renderWindow->Render();
  renderWindowInteractor->Start();
}

int main() {
  vtkNew<vtkSimplePointsReader> point_reader;
  point_reader->SetFileName(DATA_PATH "vectors.txt");
  point_reader->Update();

  mr::MeshReader mesh_reader;
  auto mesh = mesh_reader.Read(DATA_PATH "original.mesh");
  AddColor(mesh);

  mr::RefineJudger judger(point_reader->GetOutput());
  auto cell_ids_to_refine = judger.FindCellsToRefine(mesh);

  mr::InterpolateRefiner refiner(
      std::make_shared<mr::IDW>(point_reader->GetOutput()));
  auto refined_mesh = refiner.Refine(mesh, cell_ids_to_refine);
  AddColor(refined_mesh);
  /*
  {
    mr::MeshWriter writer;
    writer.Write("out.mesh", refined_mesh);
  }

  {
    mr::GRDWriter writer;
    writer.Write("out.grd", refined_mesh);
  }
*/
  Render(mesh, refined_mesh);

  return EXIT_SUCCESS;
}