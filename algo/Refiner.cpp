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

vtkNew<vtkUnsignedCharArray> Refiner::Refine(vtkPolyData* mesh) const {
  mesh->EditableOn();

  auto point_ids_to_refine = FindPointsToRefine(mesh);
  auto cell_ids_to_refine = FindUniqueCellsByPoints(mesh, point_ids_to_refine);

  auto refined_cell_ids = cell_ids_to_refine;
  std::unordered_map<vtkIdType, std::array<float, 3>> refined_boundary_cell_ids;

  std::map<std::pair<vtkIdType, vtkIdType>, vtkIdType> mid_point_id_map;

  for (auto cell_id : cell_ids_to_refine) {
    GenerateSubdivision(mesh, cell_id, refined_cell_ids,
                        refined_boundary_cell_ids, mid_point_id_map);
  }

  vtkNew<vtkUnsignedCharArray> cell_colors;
  cell_colors->SetNumberOfComponents(3);
  cell_colors->SetNumberOfTuples(mesh->GetNumberOfCells());

  for (size_t i = 0; i < mesh->GetNumberOfCells(); ++i) {
    float rgb[3] = {100, 0, 0};

    if (refined_cell_ids.find(i) != refined_cell_ids.end()) {
      rgb[1] = 100;
    } else if (auto it = refined_boundary_cell_ids.find(i);
               it != refined_boundary_cell_ids.end()) {
      rgb[0] = it->second[0];
      rgb[1] = it->second[1];
      rgb[2] = it->second[2];
    }

    cell_colors->InsertTuple(i, rgb);
  }

  mesh->GetCellData()->SetScalars(cell_colors);

  return cell_colors;
}

std::unordered_set<vtkIdType> Refiner::FindUniqueCellsByPoints(
    vtkPolyData* mesh, const std::vector<vtkIdType>& point_ids) {
  std::unordered_set<vtkIdType> res;

  for (vtkIdType point_id : point_ids) {
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

std::array<vtkIdType, 4> Refiner::GenerateSubdivision(
    vtkPolyData* mesh, vtkIdType cell_id,
    std::unordered_set<vtkIdType>& refined_cell_ids,
    std::unordered_map<vtkIdType, std::array<float, 3>>&
        refined_boundary_cell_ids,
    std::map<std::pair<vtkIdType, vtkIdType>, vtkIdType> mid_point_id_map) {
  std::array<vtkIdType, 4> res;

  const vtkIdType* point_ids;
  mesh->GetCell(cell_id, point_ids);

  double p[3][3];
  for (size_t i = 0; i < 3; ++i) mesh->GetPoint(point_ids[i], p[i]);

  vtkIdType interpolated_point_ids[3];
  for (size_t i = 0; i < 3; ++i) {
    double interpolated_point[3];
    for (size_t j = 0; j < 3; ++j) {
      interpolated_point[j] = (p[i][j] + p[(i + 1) % 3][j]) / 2;
    }

    auto two_point_ids = std::make_pair(point_ids[i], point_ids[(i + 1) % 3]);
    if (auto it = mid_point_id_map.find(two_point_ids);
        it == mid_point_id_map.end()) {
      interpolated_point_ids[i] = mid_point_id_map[two_point_ids] =
          mesh->GetPoints()->InsertNextPoint(interpolated_point);
    } else {
      interpolated_point_ids[i] = it->second;
    }
  }

  // Neighbors
  for (size_t i = 0; i < 3; ++i) {
    vtkNew<vtkIdList> cell_ids;
    mesh->GetCellEdgeNeighbors(cell_id, point_ids[i], point_ids[(i + 1) % 3],
                               cell_ids);
    for (size_t j = 0; j < cell_ids->GetNumberOfIds(); ++j) {
      auto neighbor_cell_id = cell_ids->GetId(j);
      if (refined_cell_ids.find(neighbor_cell_id) != refined_cell_ids.end())
        continue;
      if (auto it = refined_boundary_cell_ids.find(neighbor_cell_id);
          it != refined_boundary_cell_ids.end()) {
        it->second[0] = 100;
        continue;
      }

      vtkIdType n_points;
      const vtkIdType* neighbor_point_ids;
      mesh->GetCellPoints(neighbor_cell_id, n_points, neighbor_point_ids);

      assert(n_points == 3);

      vtkIdType new_point_ids[3] = {interpolated_point_ids[i]};
      for (size_t k = 0; k < 3; ++k) {
        if (neighbor_point_ids[k] != point_ids[i] &&
            neighbor_point_ids[k] != point_ids[(i + 1) % 3]) {
          new_point_ids[1] = neighbor_point_ids[k];
          break;
        }
      }

      new_point_ids[2] = point_ids[(i + 1) % 3];
      mesh->GetStrips()->ReplaceCellAtId(neighbor_cell_id, 3, new_point_ids);
      refined_boundary_cell_ids.insert({neighbor_cell_id, {0, 0, 255}});

      new_point_ids[2] = point_ids[i];
      refined_boundary_cell_ids.insert(
          {mesh->GetStrips()->InsertNextCell(3, new_point_ids), {0, 255, 255}});
    }
  }

  // Self
  mesh->GetStrips()->ReplaceCellAtId(cell_id, 3, interpolated_point_ids);
  res[0] = cell_id;

  for (size_t i = 0; i < 3; ++i) {
    refined_cell_ids.insert(res[i + 1] = mesh->GetStrips()->InsertNextCell(3));
    mesh->GetStrips()->InsertCellPoint(interpolated_point_ids[i]);
    mesh->GetStrips()->InsertCellPoint(interpolated_point_ids[(i + 1) % 3]);
    mesh->GetStrips()->InsertCellPoint(point_ids[(i + 1) % 3]);
  }

  return res;
}

}  // namespace mr
