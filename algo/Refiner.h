#pragma once

#include <vtkEdgeTable.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkStaticPointLocator2D.h>

namespace mr {

struct RefinerConfig {
  //! @brief 细分次数（每次都会将需要细分的三角形在各边中点处细分）
  int refine_times = 2;

  //! @brief 判断是否需要细分的半径（以当前判断的点为中心）
  double judge_radius = 20.0;

  //! @brief 细分后的三角形最小面积
  //!
  //! 面积小于该面积的三角形将不再细分
  double min_triangle_area = 50.0;

  //! @brief z值差阈值
  //!
  //! 如果存在一个judge_radius半径内的点，
  //! 其z值与判断中心点的差值超过该阈值，
  //! 则需要细分与该判断中心点相关联的三角形
  double delta_z_threshold = 2;

  //! @brief 为插值点计算z值时，从点数据集中采样的半径（以插值点为中心）
  double sample_radius = 20.0;

  //! @brief IDW算法的p指数参数（通过点数据集为插值点计算z值时，采用IDW算法）
  int idw_p = 2;

  //! @brief 网格线性插值计算z值的比重
  //!
  //! 为插值点计算z值时，会将通过点数据集计算的IDW插值结果与其在网格内线性插值计算的结果进行平均。
  //! 该参数决定了其在网格内线性插值计算的结果的比重。
  double original_z_weight = 0.;
};

class Refiner {
 public:
  Refiner(vtkPolyData* high_res_points, RefinerConfig config = {});

  vtkNew<vtkPolyData> Refine(vtkPolyData* mesh) const;

 private:
  RefinerConfig config_;
  vtkPolyData* high_res_points_;
  vtkNew<vtkStaticPointLocator2D> locator_;

  //! @brief
  //! @param mesh
  //! @param cell_id
  //! @param p1
  //! @param p2
  //! @param edge_data
  //! @param cell_ids To avoid frequent memory allocation
  //! @return
  static vtkIdType FindEdge(vtkPolyData* mesh, vtkIdType cell_id, vtkIdType p1,
                            vtkIdType p2, vtkIntArray* edge_data,
                            vtkIdList* cell_ids);

  vtkIdType InterpolatePosition(vtkPoints* input_points,
                                vtkPoints* output_points, vtkIdList* stencil,
                                double* weights) const;

  double GetInterpolatedZFromHighResData(double* x) const;

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

  vtkNew<vtkIdList> FindUniqueCellsByPoints(vtkPolyData* mesh,
                                            vtkIdList* point_ids) const;

  vtkNew<vtkIdList> FindPointsToRefine(vtkPolyData* mesh) const;
};

}  // namespace mr