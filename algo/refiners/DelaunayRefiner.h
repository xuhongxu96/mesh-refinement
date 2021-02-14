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

  vtkNew<vtkPolyData> Refine(vtkPolyData* input, vtkIdList* cell_ids_to_refine,
                             vtkPoints* degen_points) const override;

 private:
  std::shared_ptr<IInterpolater> interpolater_;
  DelaunayRefinerConfig config_;

  //! @brief
  //! @param input
  //! @param output
  //! @param cell_ids_to_refine
  //! @return map point id to cell id that contains it
  std::unordered_map<vtkIdType, vtkIdType> GeneratePoints(
      vtkPolyData* input, vtkPoints* output,
      vtkIdList* cell_ids_to_refine) const;
};

}  // namespace mr
