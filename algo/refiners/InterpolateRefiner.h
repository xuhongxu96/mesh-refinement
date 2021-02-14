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
  //! @brief 细分次数（每次都会将需要细分的三角形在各边中点处细分）
  int refine_times = 2;

  //! @brief 网格线性插值计算z值的比重
  //!
  //! 为插值点计算z值时，会将通过点数据集计算的IDW插值结果与其在网格内线性插值计算的结果进行平均。
  //! 该参数决定了其在网格内线性插值计算的结果的比重。
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