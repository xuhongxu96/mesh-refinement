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

static std::unordered_set<vtkIdType> highlight_cell_ids;
static std::unordered_set<vtkIdType> highlight_point_ids;

void AddColor(vtkPolyData* res,
              const std::unordered_map<vtkIdType, std::array<float, 3>>&
                  cell_color_map = {},
              const std::unordered_map<vtkIdType, std::array<float, 3>>&
                  point_color_map = {}) {
  vtkNew<vtkUnsignedCharArray> colors;
  colors->SetNumberOfComponents(3);
  colors->SetNumberOfTuples(res->GetNumberOfCells());

  for (size_t i = 0; i < res->GetNumberOfCells(); ++i) {
    bool is_point_color_mapped = false;

    vtkIdType n_pts;
    const vtkIdType* pts;
    res->GetCellPoints(i, n_pts, pts);
    for (vtkIdType j = 0; j < n_pts; ++j) {
      if (highlight_point_ids.find(pts[j]) != highlight_point_ids.end()) {
        is_point_color_mapped = true;
        colors->InsertTuple(i, std::array<float, 3>{155, 30, 255}.data());
        break;
      } else if (auto it = point_color_map.find(pts[j]);
                 it != point_color_map.end()) {
        is_point_color_mapped = true;
        colors->InsertTuple(i, it->second.data());
        break;
      }
    }

    if (is_point_color_mapped) {
      continue;
    }

    if (highlight_cell_ids.find(i) != highlight_cell_ids.end()) {
      colors->InsertTuple(i, std::array<float, 3>{155, 30, 255}.data());
    } else if (auto it = cell_color_map.find(i); it != cell_color_map.end()) {
      colors->InsertTuple(i, it->second.data());
    } else {
      float rgb[3] = {0, 0, 100};
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

vtkSmartPointer<vtkPolyData> OptimizeAngle(vtkPolyData* input, int times) {
  static mr::AngleOptimizer angle_optimizer;

  std::unordered_set<vtkIdType> optimized_cell_ids;
  std::unordered_set<vtkIdType> failed_to_optimize_cell_ids;

  vtkSmartPointer<vtkPolyData> res = input;
  for (int angle_optimize_times = 0; angle_optimize_times < times;
       ++angle_optimize_times) {
    failed_to_optimize_cell_ids.clear();
    res = angle_optimizer.Optimize(res, optimized_cell_ids,
                                   failed_to_optimize_cell_ids);
  }

  {
    std::unordered_map<vtkIdType, std::array<float, 3>> color_map;
    for (auto it : optimized_cell_ids) {
      color_map[it] = {0, 100, 0};
    }

    for (auto it : failed_to_optimize_cell_ids) {
      color_map[it] = {100, 0, 0};
    }

    AddColor(input, color_map);
    AddColor(res, color_map);
  }

  return res;
}

vtkSmartPointer<vtkPolyData> OptimizeConnection(vtkPolyData* input) {
  std::unordered_set<vtkIdType> optimized_point_ids;

  static mr::ConnectionOptimizer connection_optimizer;

  vtkSmartPointer<vtkPolyData> res =
      connection_optimizer.Optimize(input, optimized_point_ids);

  {
    std::unordered_map<vtkIdType, std::array<float, 3>> pt_color_map;
    for (auto it : optimized_point_ids) {
      pt_color_map[it] = {100, 0, 0};
    }

    AddColor(input, {}, pt_color_map);
    AddColor(res, {}, pt_color_map);
  }

  return res;
}

static auto point_reader = ([] {
  vtkNew<vtkSimplePointsReader> res;

  res->SetFileName(DATA_PATH "vectors.txt");
  res->Update();
  return res;
})();

vtkSmartPointer<vtkPolyData> Refine(vtkPolyData* mesh,
                                    vtkNew<vtkPoints>& degen_points) {
  static auto refiner = ([&] {
    return mr::DelaunayRefiner(
        std::make_shared<mr::IDW>(point_reader->GetOutput()));
  })();

  static mr::RefineJudger judger(point_reader->GetOutput());
  auto cell_ids_to_refine = judger.FindCellsToRefine(mesh);

  return refiner.Refine(mesh, cell_ids_to_refine, degen_points);
}

int main() {
  mr::MeshReader mesh_reader;
  auto mesh = mesh_reader.Read(DATA_PATH "original.mesh");
  AddColor(mesh);

  vtkNew<vtkPoints> degen_points;

  // Refine
  auto refined_mesh = Refine(mesh, degen_points);

  auto res1 = OptimizeAngle(refined_mesh, 2);
  auto res2 = OptimizeConnection(res1);
  for (int i = 0; i < 1; ++i) {
    res1 = OptimizeAngle(res2, 1);
    res2 = OptimizeConnection(res1);
  }
  refined_mesh = OptimizeAngle(res2, 1);

  res2 = Refine(res1, degen_points);

  for (int i = 0; i < 4; ++i) {
    res1 = OptimizeAngle(res2, 1);
    res2 = OptimizeConnection(res1);
  }
  res1 = OptimizeAngle(res2, 1);

  vtkPolyData* left = refined_mesh;
  vtkPolyData* right = res1;

  {
    mr::MeshWriter writer;
    writer.Write("out.mesh", right);
  }

  {
    mr::GRDWriter writer;
    writer.Write("out.grd", right);
  }

  Render(left, right, degen_points);

  return EXIT_SUCCESS;
}