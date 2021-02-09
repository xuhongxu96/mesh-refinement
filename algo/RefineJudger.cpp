#include "RefineJudger.h"

#include <vtkTriangle.h>

namespace mr {

RefineJudger::RefineJudger(vtkPolyData* high_res_points,
                           RefineJudgerConfig config)
    : high_res_points_{high_res_points}, config_{std::move(config)} {
  locator_->SetDataSet(high_res_points);
  locator_->BuildLocator();
}

vtkNew<vtkIdList> RefineJudger::FindCellsToRefine(vtkPolyData* mesh) const {
  return FindUniqueCellsByPoints(mesh, FindPointsToRefine(mesh));
}

vtkNew<vtkIdList> RefineJudger::FindUniqueCellsByPoints(
    vtkPolyData* mesh, vtkIdList* point_ids) const {
  vtkNew<vtkIdList> res;

  for (vtkIdType point_id : *point_ids) {
    vtkNew<vtkIdList> cell_ids;
    mesh->GetPointCells(point_id, cell_ids);

    for (vtkIdType j = 0; j < cell_ids->GetNumberOfIds(); ++j) {
      vtkIdType cell_id = cell_ids->GetId(j);

      vtkIdType n_points;
      const vtkIdType* cell_point_ids;
      mesh->GetCellPoints(cell_id, n_points, cell_point_ids);
      if (n_points != 3) throw std::runtime_error("Only support triangles");

      double p[3][3];
      for (int i = 0; i < 3; ++i) mesh->GetPoint(cell_point_ids[i], p[i]);
      auto area = vtkTriangle::TriangleArea(p[0], p[1], p[2]);

      if (area >= config_.min_triangle_area) {
        res->InsertUniqueId(cell_id);
      }
    }
  }

  return res;
}

vtkNew<vtkIdList> RefineJudger::FindPointsToRefine(vtkPolyData* mesh) const {
  vtkNew<vtkIdList> res;

  auto points = mesh->GetPoints();

  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
    double p[3];
    points->GetPoint(i, p);

    vtkNew<vtkIdList> point_ids_in_radius;
    locator_->FindPointsWithinRadius(config_.judge_radius, p,
                                     point_ids_in_radius);

    for (vtkIdType j = 0; j < point_ids_in_radius->GetNumberOfIds(); ++j) {
      double p_in_radius[3];
      vtkIdType point_in_radius = point_ids_in_radius->GetId(j);
      high_res_points_->GetPoint(point_in_radius, p_in_radius);

      if (abs(p[2] - p_in_radius[2]) >= config_.delta_z_threshold) {
        res->InsertNextId(i);
        break;
      }
    }
  }

  return res;
}

}  // namespace mr