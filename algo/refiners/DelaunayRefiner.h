#pragma once

#include <interpolaters/IInterpolater.h>
#include <vtkPoints.h>

#include "IRefiner.h"

namespace mr {

struct DelaunayRefinerConfig {
  double resolution_x = 10.;
  double resolution_y = 10.;
};

class DelaunayRefiner : public IRefiner {
 public:
  DelaunayRefiner(std::shared_ptr<IInterpolater> interpolater,
                  DelaunayRefinerConfig config = {});

  vtkNew<vtkPolyData> Refine(vtkPolyData* mesh,
                             vtkIdList* cell_ids_to_refine) const override;

 private:
  std::shared_ptr<IInterpolater> interpolater_;
  DelaunayRefinerConfig config_;

  void GeneratePoints(vtkPolyData* input, vtkPolyData* output,
                      vtkIdList* cell_ids_to_refine) const;
};

}  // namespace mr
