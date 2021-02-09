#pragma once

#include <vtkPolyData.h>
#include <vtkStaticPointLocator2D.h>

#include "IInterpolater.h"

namespace mr {

struct IDWConfig {
  //! @brief 为插值点计算z值时，从点数据集中采样的半径（以插值点为中心）
  double sample_radius = 20.0;

  //! @brief IDW算法的p指数参数（通过点数据集为插值点计算z值时，采用IDW算法）
  int idw_p = 2;
};

class IDW : public IInterpolater {
 public:
  IDW(vtkPolyData* sample_points, IDWConfig config = {});
  double GetValue(double x[3]) const override;

 private:
  vtkPolyData* sample_points_;
  IDWConfig config_;
  vtkNew<vtkStaticPointLocator2D> locator_;
};
}  // namespace mr