#pragma once

#include <vtkPoints.h>
#include <vtkPolyData.h>

namespace mr {

struct IRefiner {
  virtual vtkNew<vtkPolyData> Refine(vtkPolyData* mesh,
                                     vtkIdList* cell_ids_to_refine,
                                     vtkPoints* degen_points) const = 0;
};

}  // namespace mr