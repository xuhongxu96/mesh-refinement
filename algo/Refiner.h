#pragma once

#include <vtkPolyData.h>
#include <vtkStaticPointLocator2D.h>
#include <vtkUnsignedCharArray.h>

#include <array>
#include <map>
#include <unordered_map>
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

  vtkNew<vtkUnsignedCharArray> Refine(vtkPolyData* mesh) const;

 private:
  RefinerConfig config_;
  vtkPolyData* high_res_points_;
  vtkNew<vtkStaticPointLocator2D> locator_;

  static std::unordered_set<vtkIdType> FindUniqueCellsByPoints(
      vtkPolyData* mesh, const std::vector<vtkIdType>& point_ids);

  static std::array<vtkIdType, 4> GenerateSubdivision(
      vtkPolyData* mesh, vtkIdType cell_id,
      std::unordered_set<vtkIdType>& refined_cell_ids,
      std::unordered_map<vtkIdType, std::array<float, 3>>&
          refined_boundary_cell_ids,
      std::map<std::pair<vtkIdType, vtkIdType>, vtkIdType> mid_point_id_map);

  std::vector<vtkIdType> FindPointsToRefine(vtkPolyData* mesh) const;
};

}  // namespace mr