#pragma once

#include <vtkPolyData.h>

namespace mr {

struct IRefiner {
  virtual vtkNew<vtkPolyData> Refine(vtkPolyData* mesh,
                                     vtkIdList* cell_ids_to_refine) const = 0;
};

}  // namespace mr