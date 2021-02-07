#pragma once

#include <vtkEdgeTable.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkStaticPointLocator2D.h>

namespace mr {

struct RefinerConfig {
  //! @brief ϸ�ִ�����ÿ�ζ��Ὣ��Ҫϸ�ֵ��������ڸ����е㴦ϸ�֣�
  int refine_times = 2;

  //! @brief �ж��Ƿ���Ҫϸ�ֵİ뾶���Ե�ǰ�жϵĵ�Ϊ���ģ�
  double judge_radius = 20.0;

  //! @brief ϸ�ֺ����������С���
  //!
  //! ���С�ڸ�����������ν�����ϸ��
  double min_triangle_area = 50.0;

  //! @brief zֵ����ֵ
  //!
  //! �������һ��judge_radius�뾶�ڵĵ㣬
  //! ��zֵ���ж����ĵ�Ĳ�ֵ��������ֵ��
  //! ����Ҫϸ������ж����ĵ��������������
  double delta_z_threshold = 2;

  //! @brief Ϊ��ֵ�����zֵʱ���ӵ����ݼ��в����İ뾶���Բ�ֵ��Ϊ���ģ�
  double sample_radius = 20.0;

  //! @brief IDW�㷨��pָ��������ͨ�������ݼ�Ϊ��ֵ�����zֵʱ������IDW�㷨��
  int idw_p = 2;

  //! @brief �������Բ�ֵ����zֵ�ı���
  //!
  //! Ϊ��ֵ�����zֵʱ���Ὣͨ�������ݼ������IDW��ֵ������������������Բ�ֵ����Ľ������ƽ����
  //! �ò����������������������Բ�ֵ����Ľ���ı��ء�
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