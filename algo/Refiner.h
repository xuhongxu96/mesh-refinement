#pragma once

#include <vtkPolyData.h>
#include <vtkStaticPointLocator2D.h>

#include <array>
#include <unordered_set>
#include <vector>

namespace mr {

struct RefinerConfig {
  double radius = 100.0;
  double delta_z_threshold = 5;
};

class Refiner {
 public:
  Refiner(vtkPolyData* high_res_points, RefinerConfig config = {});

  void Refine(vtkPolyData* mesh) const;

 private:
  RefinerConfig config_;
  vtkPolyData* high_res_points_;
  vtkNew<vtkStaticPointLocator2D> locator_;

  static std::unordered_set<vtkIdType> FindUniqueCellsByPoints(
      vtkPolyData* mesh, const std::vector<vtkIdType>& point_ids);

  std::vector<vtkIdType> FindPointsToRefine(vtkPolyData* mesh) const;
  std::array<vtkIdType, 4> GenerateSubdivision(vtkPolyData* mesh,
                                               vtkIdType cell_id) const;
};

}  // namespace mr