#include "DelaunayRefiner.h"

#include <vtkCellArray.h>
#include <vtkDelaunay2D.h>
#include <vtkFeatureEdges.h>
#include <vtkIdFilter.h>
#include <vtkPointData.h>
#include <vtkPolygon.h>

#include <unordered_set>
namespace mr {
DelaunayRefiner::DelaunayRefiner(std::shared_ptr<IInterpolater> interpolater,
                                 DelaunayRefinerConfig config)
    : interpolater_{interpolater}, config_{std::move(config)} {}

vtkNew<vtkPolyData> DelaunayRefiner::Refine(
    vtkPolyData* mesh, vtkIdList* cell_ids_to_refine) const {
  vtkNew<vtkIdFilter> id_filter;
  id_filter->SetInputData(mesh);
  id_filter->SetPointIdsArrayName("ids");
  id_filter->SetCellIds(false);
  id_filter->Update();

  vtkNew<vtkFeatureEdges> feature_edges;
  feature_edges->SetInputConnection(id_filter->GetOutputPort());
  feature_edges->BoundaryEdgesOn();
  feature_edges->FeatureEdgesOff();
  feature_edges->ManifoldEdgesOff();
  feature_edges->NonManifoldEdgesOff();
  feature_edges->Update();

  auto boundary_arr =
      feature_edges->GetOutput()->GetPointData()->GetArray("ids");

  vtkNew<vtkPolyData> res;
  res->DeepCopy(mesh);
  res->BuildLinks();

  for (vtkIdType cell_id : *cell_ids_to_refine) {
    res->DeleteCell(cell_id);
  }

  res->RemoveDeletedCells();

  GeneratePoints(mesh, res, cell_ids_to_refine);

  vtkNew<vtkCellArray> boundary_cells;
  vtkNew<vtkPolygon> boundary_poly;
  for (vtkIdType i = 0; i < feature_edges->GetOutput()->GetNumberOfPoints();
       ++i) {
    boundary_poly->GetPointIds()->InsertNextId(i);
  }
  boundary_cells->InsertNextCell(boundary_poly);

  vtkNew<vtkPolyData> boundary;
  boundary->SetPoints(res->GetPoints());
  boundary->SetPolys(boundary_cells);

  vtkNew<vtkDelaunay2D> delaunay;
  delaunay->SetInputData(res);
  // delaunay->SetSourceData(boundary);
  delaunay->Update();

  res->ShallowCopy(delaunay->GetOutput());

  return res;
}

void DelaunayRefiner::GeneratePoints(vtkPolyData* input, vtkPolyData* output,
                                     vtkIdList* cell_ids_to_refine) const {
  for (vtkIdType cell_id : *cell_ids_to_refine) {
    double bounds[6];
    input->GetCellBounds(cell_id, bounds);

    double p[3];
    for (p[0] = bounds[0] + config_.resolution_x; p[0] < bounds[1];
         p[0] += config_.resolution_x) {
      for (p[1] = bounds[2] + config_.resolution_y; p[1] < bounds[3];
           p[1] += config_.resolution_y) {
        double z = interpolater_->GetValue(p);
        if (z == 0) continue;

        output->GetPoints()->InsertNextPoint(p[0], p[1], z);
      }
    }
  }
}

}  // namespace mr