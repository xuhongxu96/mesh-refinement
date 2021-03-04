#pragma once

#include <vtkPolyData.h>
#include <vtkStaticPointLocator2D.h>

namespace mr {

struct RefineJudgerConfig {
  //! @brief 判断是否需要细分的半径（以当前判断的点为中心）
  double judge_radius = 20.0;

  //! @brief z值差阈值
  //!
  //! 如果存在一个judge_radius半径内的点，
  //! 其z值与判断中心点的差值超过该阈值，
  //! 则需要细分与该判断中心点相关联的三角形
  double delta_z_threshold = 2;

  //! @brief 细分后的三角形最小面积
  //!
  //! 面积小于该面积的三角形将不再细分
  double min_triangle_area = 50.0;

  //! @brief 是否细分边缘
  bool refine_edge = false;
};

class RefineJudger {
 public:
  RefineJudger(vtkPolyData* high_res_points, RefineJudgerConfig config = {});

  vtkNew<vtkIdList> FindCellsToRefine(vtkPolyData* mesh) const;

 private:
  vtkPolyData* high_res_points_;
  RefineJudgerConfig config_;
  vtkNew<vtkStaticPointLocator2D> locator_;

  vtkNew<vtkIdList> FindUniqueCellsByPoints(vtkPolyData* mesh,
                                            vtkIdList* point_ids) const;

  vtkNew<vtkIdList> FindPointsToRefine(vtkPolyData* mesh) const;
};
}  // namespace mr