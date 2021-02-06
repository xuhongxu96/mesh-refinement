#pragma once

#include <vtkEdgeTable.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkStaticPointLocator2D.h>

namespace mr {

struct RefinerConfig {
  int refine_times = 1;
  double radius = 100.0;
  double delta_z_threshold = 5;
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

  static vtkIdType InterpolatePosition(vtkPoints* input_points,
                                       vtkPoints* output_points,
                                       vtkIdList* stencil, double* weights);

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

  static vtkNew<vtkIdList> FindUniqueCellsByPoints(vtkPolyData* mesh,
                                                   vtkIdList* point_ids);

  vtkNew<vtkIdList> FindPointsToRefine(vtkPolyData* mesh) const;
};

}  // namespace mr