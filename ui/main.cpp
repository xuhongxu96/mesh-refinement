#include <interpolaters/IDW.h>
#include <io/GRDWriter.h>
#include <io/MeshReader.h>
#include <io/MeshWriter.h>
#include <judgers/RefineJudger.h>
#include <optimizers/AngleOptimizer.h>
#include <refiners/DelaunayRefiner.h>
#include <refiners/InterpolateRefiner.h>
#include <vtkActor.h>
#include <vtkCellData.h>
#include <vtkDataSetMapper.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSimplePointsReader.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertexGlyphFilter.h>

#include "MouseInteractorStyle.h"

vtkNew<vtkPoints> degen_points;
std::unordered_set<vtkIdType> optimized_cell_ids;
std::unordered_set<vtkIdType> failed_to_optimize_cell_ids;

void AddColor(bool for_refined, vtkPolyData* res) {
  vtkNew<vtkUnsignedCharArray> colors;
  colors->SetNumberOfComponents(3);
  colors->SetNumberOfTuples(res->GetNumberOfCells());

  for (size_t i = 0; i < res->GetNumberOfCells(); ++i) {
    if (!for_refined) {
      float rgb[3] = {100, 0, 0};
      colors->InsertTuple(i, rgb);
      continue;
    }

    if (optimized_cell_ids.find(i) != optimized_cell_ids.end()) {
      float rgb[3] = {0, 100, 0};
      colors->InsertTuple(i, rgb);
    } else if (failed_to_optimize_cell_ids.find(i) !=
               failed_to_optimize_cell_ids.end()) {
      float rgb[3] = {0, 0, 100};
      colors->InsertTuple(i, rgb);
    } else {
      float rgb[3] = {100, 0, 0};
      colors->InsertTuple(i, rgb);
    }
  }

  res->GetCellData()->SetScalars(colors);
}

void Render(vtkPolyData* left_mesh, vtkPolyData* right_mesh,
            vtkPoints* degen_points) {
  vtkPolyData* meshes[2] = {left_mesh, right_mesh};
  double viewports[2][4] = {{0.0, 0.0, 0.5, 1.0}, {0.5, 0.0, 1.0, 1.0}};

  vtkNew<vtkRenderWindow> renderWindow;
  vtkNew<MouseInteractorStyle> style;

  for (int i = 0; i < 2; ++i) {
    vtkNew<vtkPolyData> point_mesh;
    point_mesh->SetPoints(degen_points);

    vtkNew<vtkVertexGlyphFilter> vertex_filter;
    vertex_filter->SetInputData(point_mesh);
    vertex_filter->Update();

    vtkNew<vtkPolyData> point_mesh2;
    point_mesh2->ShallowCopy(vertex_filter->GetOutput());

    vtkNew<vtkActor> point_actor;
    {
      vtkNew<vtkPolyDataMapper> mapper;
      mapper->SetInputData(point_mesh2);
      point_actor->SetMapper(mapper);
      point_actor->GetProperty()->SetPointSize(3);
      point_actor->GetProperty()->SetColor(0, 0, 0);
    }

    vtkNew<vtkActor> actor;
    {
      vtkNew<vtkPolyDataMapper> mapper;
      mapper->SetInputData(meshes[i]);
      actor->SetMapper(mapper);
      actor->GetProperty()->SetSpecular(.4);
      actor->GetProperty()->SetSpecularPower(30.0);
    }

    vtkNew<vtkRenderer> renderer;
    renderer->SetViewport(viewports[i]);

    renderWindow->AddRenderer(renderer);

    style->SetCurrentRenderer(renderer);
    style->DataMap.insert({renderer, meshes[i]});

    renderer->AddActor(actor);
    renderer->AddActor(point_actor);
    renderer->SetBackground(255, 255, 255);
  }

  renderWindow->SetWindowName("Mesh Refinement");

  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->SetInteractorStyle(style);

  renderWindow->Render();
  renderWindowInteractor->Start();
}

int main() {
  vtkNew<vtkSimplePointsReader> point_reader;
  point_reader->SetFileName(DATA_PATH "vectors.txt");
  point_reader->Update();

  mr::MeshReader mesh_reader;
  auto mesh = mesh_reader.Read(DATA_PATH "original.mesh");
  AddColor(false, mesh);

  mr::RefineJudger judger(point_reader->GetOutput());
  auto cell_ids_to_refine = judger.FindCellsToRefine(mesh);

  mr::DelaunayRefiner refiner(
      std::make_shared<mr::IDW>(point_reader->GetOutput()));

  mr::AngleOptimizer optimizer;

  auto refined_mesh =
      optimizer.Optimize(refiner.Refine(mesh, cell_ids_to_refine, degen_points),
                         optimized_cell_ids, failed_to_optimize_cell_ids);

  AddColor(true, refined_mesh);
  {
    mr::MeshWriter writer;
    writer.Write("out.mesh", refined_mesh);
  }

  {
    mr::GRDWriter writer;
    writer.Write("out.grd", refined_mesh);
  }
  Render(mesh, refined_mesh, degen_points);

  return EXIT_SUCCESS;
}