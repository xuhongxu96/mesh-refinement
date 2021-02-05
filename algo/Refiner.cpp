#include "Refiner.h"

#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkPointInterpolator.h>

namespace mr {

Refiner::Refiner(vtkPolyData* high_res_points, RefinerConfig config)
    : high_res_points_{high_res_points}, config_{config} {
  locator_->SetDataSet(high_res_points);
  locator_->BuildLocator();
}

void Refiner::Refine(vtkPolyData* mesh) const {
  auto point_ids_to_refine = FindPointsToRefine(mesh);
  auto cell_ids_to_refine = FindUniqueCellsByPoints(mesh, point_ids_to_refine);

  while (!cell_ids_to_refine.empty()) {
    auto cell_id = *cell_ids_to_refine.begin();

    GenerateSubdivision(mesh, cell_id);

    cell_ids_to_refine.erase(cell_id);
  }
}

std::unordered_set<vtkIdType> Refiner::FindUniqueCellsByPoints(
    vtkPolyData* mesh, const std::vector<vtkIdType>& point_ids) {
  std::unordered_set<vtkIdType> res;

  for (auto point_id : point_ids) {
    vtkNew<vtkIdList> cell_ids;
    mesh->GetPointCells(point_id, cell_ids);

    for (vtkIdType j = 0; j < cell_ids->GetNumberOfIds(); ++j) {
      res.insert(cell_ids->GetId(j));
    }
  }

  return res;
}

std::vector<vtkIdType> Refiner::FindPointsToRefine(vtkPolyData* mesh) const {
  std::vector<vtkIdType> res;

  auto points = mesh->GetPoints();

  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
    double p[3];
    points->GetPoint(i, p);

    vtkNew<vtkIdList> point_ids_in_radius;
    locator_->FindPointsWithinRadius(config_.radius, p, point_ids_in_radius);

    for (vtkIdType j = 0; j < point_ids_in_radius->GetNumberOfIds(); ++j) {
      double p_in_radius[3];
      vtkIdType point_in_radius = point_ids_in_radius->GetId(j);
      high_res_points_->GetPoint(point_in_radius, p_in_radius);

      if (abs(p[2] - p_in_radius[2]) >= config_.delta_z_threshold) {
        res.push_back(i);
        break;
      }
    }
  }

  return res;
}

std::array<vtkIdType, 4> Refiner::GenerateSubdivision(vtkPolyData* mesh,
                                                      vtkIdType cell_id) const {
  const vtkIdType* point_ids;
  mesh->GetCell(cell_id, point_ids);

  double p[3][3];
  for (auto i = 0; i < 3; ++i) mesh->GetPoint(point_ids[i], p[i]);

  double interpolated_p[3][3];
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 3; ++j) {
      interpolated_p[i][j] = (p[i][j] + p[(i + 1) % 3][j]) / 2;
    }

    mesh->GetPoints()->InsertNextPoint(interpolated_p[i]);
  }
}

}  // namespace mr
