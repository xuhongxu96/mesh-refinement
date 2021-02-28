#include "ConnectionOptimizer.h"

#include <unordered_set>

namespace mr {

ConnectionOptimizer::ConnectionOptimizer(ConnectionOptimizerConfig config)
    : config_{std::move(config)} {}

vtkNew<vtkIdList> ConnectionOptimizer::FindPointIdsToOptimize(
    vtkPolyData* input) const {
  vtkNew<vtkIdList> res;

  vtkIdType n_points = input->GetNumberOfPoints();
  for (vtkIdType i = 0; i < n_points; ++i) {
    vtkIdType n_cells;
    vtkIdType* cell_ids;
    input->GetPointCells(i, n_cells, cell_ids);

    if (n_cells > config_.max_connection) {
      res->InsertNextId(i);
    }
  }

  return res;
}

vtkNew<vtkPolyData> ConnectionOptimizer::Optimize(
    vtkPolyData* input, std::unordered_set<vtkIdType>& optimized_cell_ids,
    std::unordered_set<vtkIdType>& failed_point_ids) const {
  vtkNew<vtkPolyData> res;
  res->EditableOn();
  res->DeepCopy(input);
  res->BuildLinks();

  vtkNew<vtkIdList> point_ids_to_optimize = FindPointIdsToOptimize(input);
  for (auto id : *point_ids_to_optimize) {
    failed_point_ids.insert(id);
  }
  // TODO: optimize connection

  // For n times, n == actual_connection - max_connection,
  // Find longest edge and half-split it,
  // make neighbor 2 points disconnect with mid point but connect to the half
  // point, make neighbor-of-neighbor 2 points connect with the half point.

  res->BuildLinks();
  res->Squeeze();
  return res;
}

}  // namespace mr