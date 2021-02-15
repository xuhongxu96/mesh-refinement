#pragma once

#include <vtkIdList.h>
#include <vtkPolyData.h>

#include <unordered_set>

namespace mr {

struct AngleOptimizerConfig {
  double max_angle = 130.;
  vtkIdType max_connection = 8;
};

class AngleOptimizer {
 public:
  AngleOptimizer(AngleOptimizerConfig config = {});

  vtkNew<vtkPolyData> Optimize(
      vtkPolyData* input, std::unordered_set<vtkIdType>& optimized_cell_ids,
      std::unordered_set<vtkIdType>& failed_cell_ids) const;

 private:
  AngleOptimizerConfig config_;
};
}  // namespace mr