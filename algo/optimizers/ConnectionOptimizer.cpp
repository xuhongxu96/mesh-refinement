#include "ConnectionOptimizer.h"

#include <vtkMath.h>

#include <array>
#include <optional>

namespace mr {
struct PointInfo {
  vtkIdType id = -1;
  double x[3];
};

struct NeighborInfo {
  PointInfo pivot;

  PointInfo farest_point;

  PointInfo l2_points[2];
  vtkIdType l2_cell_ids[2];

  PointInfo l3_points[2];
  vtkIdType l3_cell_ids[2];
};

ConnectionOptimizer::ConnectionOptimizer(ConnectionOptimizerConfig config)
    : config_{std::move(config)} {}

std::unordered_map<vtkIdType, vtkIdType>
ConnectionOptimizer::FindPointIdsWithConnectionToOptimize(
    vtkPolyData* input) const {
  std::unordered_map<vtkIdType, vtkIdType> res;

  vtkIdType n_points = input->GetNumberOfPoints();
  for (vtkIdType i = 0; i < n_points; ++i) {
    vtkIdType n_cells;
    vtkIdType* cell_ids;
    input->GetPointCells(i, n_cells, cell_ids);

    if (n_cells > config_.max_connection) {
      res[i] = n_cells;
    }
  }

  return res;
}

using id_point_map_t = std::unordered_map<vtkIdType, std::array<double, 3>>;
static id_point_map_t FindNeighborPoints(vtkPolyData* mesh,
                                         vtkIdType point_id) {
  vtkIdType n_cells;
  vtkIdType* cell_ids;
  mesh->GetPointCells(point_id, n_cells, cell_ids);

  id_point_map_t res;
  for (vtkIdType i = 0; i < n_cells; ++i) {
    vtkIdType n_points;
    const vtkIdType* points;
    mesh->GetCellPoints(cell_ids[i], n_points, points);

    for (vtkIdType j = 0; j < n_points; ++j) {
      mesh->GetPoint(points[j], res[points[j]].data());
    }
  }

  return res;
}

static bool GetNeighborInfo(vtkPolyData* input, vtkIdType point_id,
                            vtkIdType nei_point_id, NeighborInfo& res) {
  {
    vtkNew<vtkIdList> l2_neighbor_cell_ids;
    input->GetCellEdgeNeighbors(-1, point_id, nei_point_id,
                                l2_neighbor_cell_ids);

    if (l2_neighbor_cell_ids->GetNumberOfIds() != 2) return false;

    for (int i = 0; i < 2; ++i)
      res.l2_cell_ids[i] = l2_neighbor_cell_ids->GetId(i);

    for (vtkIdType i = 0; i < 2; ++i) {
      vtkIdType n_pts;
      const vtkIdType* pts;
      input->GetCellPoints(l2_neighbor_cell_ids->GetId(i), n_pts, pts);
      for (vtkIdType j = 0; j < n_pts; ++j) {
        if (pts[j] != res.pivot.id && pts[j] != nei_point_id) {
          res.l2_points[i].id = pts[j];
          input->GetPoint(pts[j], res.l2_points[i].x);
          break;
        }
      }
    }
  }

  {
    for (vtkIdType i = 0; i < 2; ++i) {
      vtkNew<vtkIdList> l3_neighbor_cell_ids;
      input->GetCellEdgeNeighbors(res.l2_cell_ids[i], point_id,
                                  res.l2_points[i].id, l3_neighbor_cell_ids);
      if (l3_neighbor_cell_ids->GetNumberOfIds() != 1) return false;

      res.l3_cell_ids[i] = l3_neighbor_cell_ids->GetId(0);

      vtkIdType n_pts;
      const vtkIdType* pts;
      input->GetCellPoints(res.l3_cell_ids[i], n_pts, pts);
      for (vtkIdType j = 0; j < n_pts; ++j) {
        if (pts[j] != res.pivot.id && pts[j] != res.l2_points[i].id) {
          res.l3_points[i].id = pts[j];
          input->GetPoint(pts[j], res.l3_points[i].x);
          break;
        }
      }
    }
  }

  return true;
}

static NeighborInfo GetNeighborInfo(vtkPolyData* input, vtkIdType point_id) {
  NeighborInfo res;

  res.pivot.id = point_id;
  input->GetPoint(point_id, res.pivot.x);

  auto neighbor_points = FindNeighborPoints(input, point_id);

  // Find longest edge
  double max_distance = 0;
  for (auto& [neighbor_point_id, neighbor_point] : neighbor_points) {
    auto distance =
        vtkMath::Distance2BetweenPoints(res.pivot.x, neighbor_point.data());
    NeighborInfo temp = res;
    if (max_distance <= distance &&
        GetNeighborInfo(input, point_id, neighbor_point_id, temp)) {
      res = temp;
      max_distance = distance;
      res.farest_point.id = neighbor_point_id;
      for (int i = 0; i < 3; ++i) res.farest_point.x[i] = neighbor_point[i];
    }
  }

  if (res.farest_point.id == -1 || max_distance == 0) {
    throw std::runtime_error("Failed to find longest edge (shouldn't happen)");
  }

  return res;
}

vtkNew<vtkPolyData> ConnectionOptimizer::Optimize(
    vtkPolyData* input,
    std::unordered_set<vtkIdType>& optimized_point_ids) const {
  optimized_point_ids.clear();

  vtkNew<vtkPolyData> res;
  res->EditableOn();
  res->DeepCopy(input);
  res->BuildLinks();

  auto point_ids_with_connection_to_optimize =
      FindPointIdsWithConnectionToOptimize(res);
  for (auto [point_id, connection] : point_ids_with_connection_to_optimize) {
    auto delta = connection - config_.max_connection;

    // For delta times
    for (int times = 0; times < delta; ++times) {
      auto info = GetNeighborInfo(res, point_id);

      // half-split longest edge
      double mid_point[3];
      for (int i = 0; i < 3; ++i)
        mid_point[i] = (info.farest_point.x[i] + info.pivot.x[i]) / 2.;

      // auto mid_point_id = res->GetPoints()->InsertNextPoint(mid_point);
      auto mid_point_id = res->InsertNextLinkedPoint(mid_point, 0);

      for (int i = 0; i < 2; ++i) {
        // replace l2 cell (pivot -> mid)
        res->RemoveReferenceToCell(info.pivot.id, info.l2_cell_ids[i]);
        res->ReplaceCellPoint(info.l2_cell_ids[i], info.pivot.id, mid_point_id);
        res->ResizeCellList(mid_point_id, 1);
        res->AddReferenceToCell(mid_point_id, info.l2_cell_ids[i]);

        // replace l3 cell (pivot -> mid)
        res->RemoveReferenceToCell(info.pivot.id, info.l3_cell_ids[i]);
        res->ReplaceCellPoint(info.l3_cell_ids[i], info.pivot.id, mid_point_id);
        res->ResizeCellList(mid_point_id, 1);
        res->AddReferenceToCell(mid_point_id, info.l3_cell_ids[i]);

        // new cell (mid, pivot, l3)
        vtkIdType pts[3] = {mid_point_id, info.pivot.id, info.l3_points[i].id};
        auto new_cell_id = res->InsertNextLinkedCell(VTK_TRIANGLE, 3, pts);
      }
    }

    optimized_point_ids.insert(point_id);
  }

  res->Squeeze();
  return res;
}

}  // namespace mr