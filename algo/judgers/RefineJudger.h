#pragma once

#include <vtkPolyData.h>
#include <vtkStaticPointLocator2D.h>

namespace mr {

struct RefineJudgerConfig {
  //! @brief �ж��Ƿ���Ҫϸ�ֵİ뾶���Ե�ǰ�жϵĵ�Ϊ���ģ�
  double judge_radius = 20.0;

  //! @brief zֵ����ֵ
  //!
  //! �������һ��judge_radius�뾶�ڵĵ㣬
  //! ��zֵ���ж����ĵ�Ĳ�ֵ��������ֵ��
  //! ����Ҫϸ������ж����ĵ��������������
  double delta_z_threshold = 2;

  //! @brief ϸ�ֺ����������С���
  //!
  //! ���С�ڸ�����������ν�����ϸ��
  double min_triangle_area = 50.0;

  //! @brief �Ƿ�ϸ�ֱ�Ե
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