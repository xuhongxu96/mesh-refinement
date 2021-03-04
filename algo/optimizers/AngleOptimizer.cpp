#include "AngleOptimizer.h"

#include <vtkMath.h>

#include <iostream>
#include <optional>
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

static int GetMaxRadianIndex(double p[3][3], double radian) {
  double vec[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      vec[i][j] = p[i][j] - p[(i + 1) % 3][j];
    }
    vec[i][2] = 0;  // calculate angle in 2d only
  }

  int max_radian_index = -1;
  for (int i = 0; i < 3; ++i) {
    double current_radian =
        vtkMath::Pi() - vtkMath::AngleBetweenVectors(vec[i], vec[(i + 1) % 3]);
    if (current_radian >= radian) max_radian_index = i;
  }

  return max_radian_index;
}

static std::optional<CellEdge> DoesCellNeedOptimization(vtkPolyData* mesh,
                                                        vtkIdType cell_id,
                                                        double radian) {
  vtkIdType n_points;
  const vtkIdType* point_ids;
  mesh->GetCellPoints(cell_id, n_points, point_ids);
  assert(n_points == 3);

  double p[3][3];
  for (int i = 0; i < 3; ++i) mesh->GetPoint(point_ids[i], p[i]);

  int max_radian_index = GetMaxRadianIndex(p, radian);
  if (max_radian_index == -1) return {};

  CellEdge res;
  res.max_radian_index = max_radian_index;

  res.max_p_id = point_ids[(max_radian_index + 1) % 3];
  res.p1_id = point_ids[max_radian_index];
  res.p2_id = point_ids[(max_radian_index + 2) % 3];

  for (int i = 0; i < 3; ++i) {
    res.max_p[i] = p[(max_radian_index + 1) % 3][i];
    res.p1[i] = p[max_radian_index][i];
    res.p2[i] = p[(max_radian_index + 2) % 3][i];
  }

  return res;
}

static std::unordered_map<vtkIdType, CellEdge> FindCellEdgeWithBigAngle(
    vtkPolyData* mesh, double radian) {
  std::unordered_map<vtkIdType, CellEdge> res;

  for (vtkIdType cell_id = 0; cell_id < mesh->GetNumberOfCells(); ++cell_id) {
    auto info = DoesCellNeedOptimization(mesh, cell_id, radian);
    if (info) res[cell_id] = *info;
  }

  return res;
}

static vtkIdType GetOppositePointId(vtkPolyData* mesh, vtkIdType cell_id,
                                    vtkIdType p1, vtkIdType p2) {
  vtkIdType n_points;
  const vtkIdType* point_ids;
  mesh->GetCellPoints(cell_id, n_points, point_ids);
  assert(n_points == 3);
  for (int i = 0; i < 3; ++i) {
    if (point_ids[i] != p1 && point_ids[i] != p2) {
      return point_ids[i];
    }
  }

  throw std::runtime_error(
      "cannot get opposite point id, which shouldn't happen");
}

static void SwapDiagonal(vtkPolyData* mesh, vtkIdType cell_id,
                         vtkIdType nei_cell_id, vtkIdType nei_point_id,
                         const CellEdge& info) {
  mesh->RemoveReferenceToCell(info.p1_id, cell_id);
  mesh->RemoveReferenceToCell(info.p2_id, nei_cell_id);
  mesh->ResizeCellList(nei_point_id, 1);
  mesh->AddReferenceToCell(nei_point_id, cell_id);
  mesh->ResizeCellList(info.max_p_id, 1);
  mesh->AddReferenceToCell(info.max_p_id, nei_cell_id);

  {
    vtkIdType new_point_ids[3] = {info.p1_id, info.max_p_id, nei_point_id};
    mesh->ReplaceCell(cell_id, 3, new_point_ids);
  }

  {
    vtkIdType new_point_ids[3] = {info.max_p_id, info.p2_id, nei_point_id};
    mesh->ReplaceCell(nei_cell_id, 3, new_point_ids);
  }
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
    if (optimized_cell_ids.find(cell_id) != optimized_cell_ids.end()) continue;

    vtkNew<vtkIdList> nei_cell_ids;
    res->GetCellEdgeNeighbors(cell_id, info.p1_id, info.p2_id, nei_cell_ids);
    if (cell_id == 188788) {
      std::cout << "a" << std::endl;
    }
    switch (nei_cell_ids->GetNumberOfIds()) {
      case 0:
        // edge cell, directly split the edge
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
        // 2. split the edge and the neighbor, need to check if neighbor is
        // really optmized
        {
          vtkIdType nei_cell_id = nei_cell_ids->GetId(0);
          vtkIdType nei_point_id =
              GetOppositePointId(res, nei_cell_id, info.p1_id, info.p2_id);

          if (optimized_cell_ids.find(nei_cell_id) != optimized_cell_ids.end())
            continue;

          bool succ = true;

          {
            double p[3][3];
            for (int i = 0; i < 3; ++i) {
              p[0][i] = info.p1[i];
              p[1][i] = info.max_p[i];
              p[2][i] = res->GetPoint(nei_point_id)[i];
            }
            if (GetMaxRadianIndex(p, max_radian) != -1) succ = false;
          }

          {
            double p[3][3];
            for (int i = 0; i < 3; ++i) {
              p[0][i] = info.max_p[i];
              p[1][i] = info.p2[i];
              p[2][i] = res->GetPoint(nei_point_id)[i];
            }
            if (GetMaxRadianIndex(p, max_radian) != -1) succ = false;
          }

          if (succ) {
            SwapDiagonal(res, cell_id, nei_cell_id, nei_point_id, info);
            optimized_cell_ids.insert(cell_id);
            optimized_cell_ids.insert(nei_cell_id);
          } else {
            double mid_p[3];
            for (int i = 0; i < 3; ++i) {
              mid_p[i] = (info.p1[i] + info.p2[i]) / 2.;
            }

            vtkIdType mid_p_id = res->InsertNextLinkedPoint(mid_p, 1);

            {
              vtkIdType new_point_ids[3] = {info.max_p_id, mid_p_id,
                                            info.p1_id};
              res->RemoveReferenceToCell(info.p2_id, cell_id);
              res->ReplaceCell(cell_id, 3, new_point_ids);
              res->ResizeCellList(mid_p_id, 1);
              res->AddReferenceToCell(mid_p_id, cell_id);
              optimized_cell_ids.insert(cell_id);
            }

            {
              vtkIdType new_point_ids[3] = {info.max_p_id, info.p2_id,
                                            mid_p_id};
              vtkIdType new_cell_id =
                  res->InsertNextLinkedCell(VTK_POLYGON, 3, new_point_ids);
              optimized_cell_ids.insert(new_cell_id);
            }

            {
              vtkIdType new_point_ids[3] = {nei_point_id, mid_p_id, info.p1_id};
              res->RemoveReferenceToCell(info.p2_id, nei_cell_id);
              res->ReplaceCell(nei_cell_id, 3, new_point_ids);
              res->ResizeCellList(mid_p_id, 1);
              res->AddReferenceToCell(mid_p_id, nei_cell_id);
              optimized_cell_ids.insert(nei_cell_id);
            }

            {
              vtkIdType new_point_ids[3] = {nei_point_id, info.p2_id, mid_p_id};
              vtkIdType new_cell_id =
                  res->InsertNextLinkedCell(VTK_POLYGON, 3, new_point_ids);
              optimized_cell_ids.insert(new_cell_id);
            }
          }
        }
        break;
      default:
        // throw std::runtime_error("neighbor cannot be greater than 1");
        failed_cell_ids.insert(cell_id);
        break;
    }
  }

  res->BuildLinks();
  res->Squeeze();
  return res;
}
}  // namespace mr