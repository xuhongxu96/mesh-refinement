#pragma once

#include <vtkPolyData.h>
#include <vtkStaticPointLocator2D.h>

#include "IInterpolater.h"

namespace mr {

struct IDWConfig {
  //! @brief Ϊ��ֵ�����zֵʱ���ӵ����ݼ��в����İ뾶���Բ�ֵ��Ϊ���ģ�
  double sample_radius = 20.0;

  //! @brief IDW�㷨��pָ��������ͨ�������ݼ�Ϊ��ֵ�����zֵʱ������IDW�㷨��
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