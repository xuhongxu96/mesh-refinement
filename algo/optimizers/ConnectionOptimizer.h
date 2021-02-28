#pragma once

#include <vtkIdList.h>
#include <vtkPolyData.h>

#include <unordered_set>

namespace mr {

struct ConnectionOptimizerConfig {
  int max_connection = 8;
};

class ConnectionOptimizer {
 public:
  ConnectionOptimizer(ConnectionOptimizerConfig config = {});

  vtkNew<vtkPolyData> Optimize(
      vtkPolyData* input, std::unordered_set<vtkIdType>& optimized_cell_ids,
      std::unordered_set<vtkIdType>& failed_point_ids) const;

 private:
  ConnectionOptimizerConfig config_;

  vtkNew<vtkIdList> FindPointIdsToOptimize(vtkPolyData* input) const;
};
}  // namespace mr
