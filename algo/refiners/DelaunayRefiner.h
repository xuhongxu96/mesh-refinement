#pragma once

#include <interpolaters/IInterpolater.h>
#include <vtkPoints.h>

#include "IRefiner.h"

namespace mr {

struct DelaunayRefinerConfig {};

class DelaunayRefiner : public IRefiner {
 public:
  DelaunayRefiner(std::shared_ptr<IInterpolater> interpolater,
                  DelaunayRefinerConfig config = {});

  vtkNew<vtkPolyData> Refine(vtkPolyData* input,
                             vtkIdList* cell_ids_to_refine) const override;

 private:
  std::shared_ptr<IInterpolater> interpolater_;
  DelaunayRefinerConfig config_;

  void GeneratePoints(vtkPolyData* input, vtkPolyData* output,
                      vtkIdList* cell_ids_to_refine) const;
};

}  // namespace mr
