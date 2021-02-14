#pragma once

#include <interpolaters/IInterpolater.h>
#include <vtkEdgeTable.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>

#include <memory>

#include "IRefiner.h"

namespace mr {

struct InterpolateRefinerConfig {
  //! @brief ϸ�ִ�����ÿ�ζ��Ὣ��Ҫϸ�ֵ��������ڸ����е㴦ϸ�֣�
  int refine_times = 2;

  //! @brief �������Բ�ֵ����zֵ�ı���
  //!
  //! Ϊ��ֵ�����zֵʱ���Ὣͨ�������ݼ������IDW��ֵ������������������Բ�ֵ����Ľ������ƽ����
  //! �ò����������������������Բ�ֵ����Ľ���ı��ء�
  double original_z_weight = 0.;
};

class InterpolateRefiner : public IRefiner {
 public:
  InterpolateRefiner(std::shared_ptr<IInterpolater> interpolater,
                     InterpolateRefinerConfig config = {});

  vtkNew<vtkPolyData> Refine(vtkPolyData* mesh, vtkIdList* cell_ids_to_refine,
                             vtkPoints* degen_points) const override;

 private:
  std::shared_ptr<IInterpolater> interpolater_;
  InterpolateRefinerConfig config_;

  vtkIdType InterpolatePosition(vtkPoints* input_points,
                                vtkPoints* output_points, vtkIdList* stencil,
                                double* weights) const;

  void GenerateSubdivisionPoints(vtkPolyData* input_ds,
                                 vtkEdgeTable* edge_table,
                                 vtkIntArray* edge_data,
                                 vtkPoints* output_points,
                                 vtkPointData* output_pd,
                                 vtkIdList* cell_ids_to_refine) const;

  void GenerateSubdivisionCells(vtkPolyData* input_ds, vtkEdgeTable* edge_table,
                                vtkIntArray* edge_data,
                                vtkCellArray* output_cells,
                                vtkCellData* output_cd) const;
};

}  // namespace mr