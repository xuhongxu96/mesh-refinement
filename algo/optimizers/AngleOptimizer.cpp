#include "AngleOptimizer.h"

#include <vtkMath.h>

#include <unordered_set>

namespace mr {
AngleOptimizer::AngleOptimizer(AngleOptimizerConfig config)
    : config_{std::move(config)} {
  if (config.max_angle <= 90) {
    throw std::runtime_error("max_angle must be greater than 90 degrees");
  }
}

static vtkIdType GetPointConnection(vtkPolyData* mesh, vtkIdType point_id) {
  vtkIdType n_connected_cells;
  vtkIdType* connected_cell_ids;
  mesh->GetPointCells(point_id, n_connected_cells, connected_cell_ids);

  std::unordered_set<vtkIdType> connected_points;
  for (vtkIdType i = 0; i < n_connected_cells; ++i) {
    vtkIdType n_points;
    const vtkIdType* point_ids;
    mesh->GetCellPoints(connected_cell_ids[i], n_points, point_ids);
    assert(n_points == 3);

    for (vtkIdType j = 0; j < n_points; ++j) {
      if (point_ids[j] != point_id) {
        connected_points.insert(point_ids[j]);
      }
    }
  }

  return connected_points.size();
}

struct CellEdge {
  vtkIdType max_p_id;
  vtkIdType p1_id;
  vtkIdType p2_id;

  double max_p[3];
  double p1[3];
  double p2[3];

  int max_radian_index;
};

static std::unordered_map<vtkIdType, CellEdge> FindCellEdgeWithBigAngle(
    vtkPolyData* mesh, double radian) {
  std::unordered_map<vtkIdType, CellEdge> res;
  for (vtkIdType cell_id = 0; cell_id < mesh->GetNumberOfCells(); ++cell_id) {
    vtkIdType n_points;
    const vtkIdType* point_ids;
    mesh->GetCellPoints(cell_id, n_points, point_ids);
    assert(n_points == 3);

    double p[3][3];
    for (int i = 0; i < 3; ++i) mesh->GetPoint(point_ids[i], p[i]);

    double vec[3][3];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        vec[i][j] = p[i][j] - p[(i + 1) % 3][j];
      }
    }

    double radians[3];
    int max_radian_index = -1;
    for (int i = 0; i < 3; ++i) {
      radians[i] = vtkMath::Pi() -
                   vtkMath::AngleBetweenVectors(vec[i], vec[(i + 1) % 3]);
      if (radians[i] >= radian) max_radian_index = i;
    }

    if (max_radian_index == -1) continue;

    res[cell_id].max_radian_index = max_radian_index;

    res[cell_id].max_p_id = point_ids[(max_radian_index + 1) % 3];
    res[cell_id].p1_id = point_ids[max_radian_index];
    res[cell_id].p2_id = point_ids[(max_radian_index + 2) % 3];

    for (int i = 0; i < 3; ++i) {
      res[cell_id].max_p[i] = p[(max_radian_index + 1) % 3][i];
      res[cell_id].p1[i] = p[max_radian_index][i];
      res[cell_id].p2[i] = p[(max_radian_index + 2) % 3][i];
    }
  }

  return res;
}

vtkNew<vtkPolyData> AngleOptimizer::Optimize(
    vtkPolyData* input, std::unordered_set<vtkIdType>& optimized_cell_ids,
    std::unordered_set<vtkIdType>& failed_cell_ids) const {
  vtkNew<vtkPolyData> res;
  res->EditableOn();
  res->DeepCopy(input);
  res->BuildLinks();

  double max_radian = vtkMath::RadiansFromDegrees(config_.max_angle);
  auto cell_edges = FindCellEdgeWithBigAngle(res, max_radian);

  for (auto& [cell_id, info] : cell_edges) {
    vtkIdType connection = GetPointConnection(res, info.max_p_id);
    if (connection >= config_.max_connection) {
      failed_cell_ids.insert(cell_id);
      continue;
    }

    vtkNew<vtkIdList> nei_cell_ids;
    res->GetCellEdgeNeighbors(cell_id, info.p1_id, info.p2_id, nei_cell_ids);
    switch (nei_cell_ids->GetNumberOfIds()) {
      case 0:
        // directly split the edge
        {
          double mid_p[3];
          for (int i = 0; i < 3; ++i) {
            mid_p[i] = (info.p1[i] + info.p2[i]) / 2.;
          }

          vtkIdType mid_p_id = res->InsertNextLinkedPoint(mid_p, 1);

          {
            vtkIdType new_point_ids[3] = {info.max_p_id, mid_p_id, info.p1_id};
            res->RemoveReferenceToCell(info.p2_id, cell_id);
            res->ReplaceCell(cell_id, 3, new_point_ids);
            res->ResizeCellList(mid_p_id, 1);
            res->AddReferenceToCell(mid_p_id, cell_id);
            optimized_cell_ids.insert(cell_id);
          }

          {
            vtkIdType new_point_ids[3] = {info.max_p_id, info.p2_id, mid_p_id};
            vtkIdType new_cell_id =
                res->InsertNextLinkedCell(VTK_POLYGON, 3, new_point_ids);
            optimized_cell_ids.insert(new_cell_id);
          }
        }
        break;
      case 1:
        // two ways:
        // 1. swap diagonal, 
        // 2. split the edge and the neighbor, need to check if neighbor is really
        // optmized
        failed_cell_ids.insert(cell_id);
        break;
      default:
        throw std::runtime_error("neighbor cannot be greater than 1");
    }
  }

  res->BuildLinks();
  res->Squeeze();
  return res;
}
}  // namespace mr