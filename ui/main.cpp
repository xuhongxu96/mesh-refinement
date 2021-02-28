#include <interpolaters/IDW.h>
#include <io/GRDWriter.h>
#include <io/MeshReader.h>
#include <io/MeshWriter.h>
#include <judgers/RefineJudger.h>
#include <optimizers/AngleOptimizer.h>
#include <optimizers/ConnectionOptimizer.h>
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

#include <array>
#include <vector>

#include "MouseInteractorStyle.h"

void AddColor(vtkPolyData* res,
              const std::unordered_map<vtkIdType, std::array<float, 3>>&
                  cell_color_map = {},
              const std::unordered_map<vtkIdType, std::array<float, 3>>&
                  point_color_map = {}) {
  vtkNew<vtkUnsignedCharArray> colors;
  colors->SetNumberOfComponents(3);
  colors->SetNumberOfTuples(res->GetNumberOfCells());

  for (size_t i = 0; i < res->GetNumberOfCells(); ++i) {
    if (auto it = cell_color_map.find(i); it != cell_color_map.end()) {
      colors->InsertTuple(i, it->second.data());
    } else {
      bool is_point_color_mapped = false;

      vtkIdType n_pts;
      const vtkIdType* pts;
      res->GetCellPoints(i, n_pts, pts);
      for (vtkIdType j = 0; j < n_pts; ++j) {
        if (auto it = point_color_map.find(pts[j]);
            it != point_color_map.end()) {
          is_point_color_mapped = true;
          colors->InsertTuple(i, it->second.data());
          break;
        }
      }

      if (!is_point_color_mapped) {
        float rgb[3] = {0, 0, 100};
        colors->InsertTuple(i, rgb);
      }
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
  AddColor(mesh);

  mr::RefineJudger judger(point_reader->GetOutput());
  auto cell_ids_to_refine = judger.FindCellsToRefine(mesh);

  mr::DelaunayRefiner refiner(
      std::make_shared<mr::IDW>(point_reader->GetOutput()));

  mr::AngleOptimizer angle_optimizer;

  vtkNew<vtkPoints> degen_points;
  std::unordered_set<vtkIdType> optimized_cell_ids;
  std::unordered_set<vtkIdType> failed_to_optimize_cell_ids;

  auto refined_mesh = refiner.Refine(mesh, cell_ids_to_refine, degen_points);
  auto optimized_mesh = angle_optimizer.Optimize(
      refined_mesh, optimized_cell_ids, failed_to_optimize_cell_ids);
  failed_to_optimize_cell_ids.clear();
  auto optimized_mesh2 = angle_optimizer.Optimize(
      optimized_mesh, optimized_cell_ids, failed_to_optimize_cell_ids);

  {
    std::unordered_map<vtkIdType, std::array<float, 3>> color_map;
    for (auto it : optimized_cell_ids) {
      color_map[it] = {0, 100, 0};
    }

    for (auto it : failed_to_optimize_cell_ids) {
      color_map[it] = {100, 0, 0};
    }

    AddColor(refined_mesh, color_map);
    AddColor(optimized_mesh2, color_map);
  }

  mr::ConnectionOptimizer connection_optimizer;

  optimized_cell_ids.clear();
  std::unordered_set<vtkIdType> failed_to_optimize_point_ids;

  auto connection_optimized_mesh = connection_optimizer.Optimize(
      optimized_mesh2, optimized_cell_ids, failed_to_optimize_point_ids);

  {
    std::unordered_map<vtkIdType, std::array<float, 3>> color_map;
    for (auto it : optimized_cell_ids) {
      color_map[it] = {0, 100, 0};
    }

    std::unordered_map<vtkIdType, std::array<float, 3>> pt_color_map;
    for (auto it : failed_to_optimize_point_ids) {
      pt_color_map[it] = {100, 0, 0};
    }

    AddColor(optimized_mesh2, color_map, pt_color_map);
    AddColor(connection_optimized_mesh, color_map, pt_color_map);
  }

  {
    mr::MeshWriter writer;
    writer.Write("out.mesh", optimized_mesh2);
  }

  {
    mr::GRDWriter writer;
    writer.Write("out.grd", optimized_mesh2);
  }
  Render(optimized_mesh2, connection_optimized_mesh, degen_points);

  return EXIT_SUCCESS;
}