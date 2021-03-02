#pragma once

#include <vtkIdList.h>
#include <vtkPolyData.h>

#include <unordered_map>
#include <unordered_set>

namespace mr {

struct ConnectionOptimizerConfig {
  int max_connection = 8;
};

class ConnectionOptimizer {
 public:
  ConnectionOptimizer(ConnectionOptimizerConfig config = {});

  vtkNew<vtkPolyData> Optimize(
      vtkPolyData* input,
      std::unordered_set<vtkIdType>& optimized_point_ids) const;

 private:
  ConnectionOptimizerConfig config_;

  std::unordered_map<vtkIdType, vtkIdType> FindPointIdsWithConnectionToOptimize(
      vtkPolyData* input) const;
};
}  // namespace mr
