#include "IDW.h"

namespace mr {

static constexpr double EPSILON = 0.001;

IDW::IDW(vtkPolyData* sample_points, IDWConfig config)
    : config_{std::move(config)}, sample_points_{sample_points} {
  locator_->SetDataSet(sample_points);
  locator_->BuildLocator();
}

double IDW::GetValue(double x[3]) const {
  double expected_z = 0;
  double sum = 0.;
  vtkNew<vtkIdList> point_ids_in_radius;
  locator_->FindPointsWithinRadius(config_.sample_radius, x,
                                   point_ids_in_radius);
  if (point_ids_in_radius->GetNumberOfIds() > 0) {
    for (auto point_id_in_radius : *point_ids_in_radius) {
      double p[3];
      sample_points_->GetPoint(point_id_in_radius, p);
      double d =
          sqrt((p[0] - x[0]) * (p[0] - x[0]) + (p[1] - x[1]) * (p[1] - x[1]));

      if (d < EPSILON) continue;

      double d_p = pow(d, config_.idw_p);

      expected_z += p[2] / d_p;
      sum += 1. / d_p;
    }

    expected_z /= sum;
  }

  return expected_z >= EPSILON ? expected_z : 0;
}

}  // namespace mr