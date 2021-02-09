#pragma once

namespace mr {
struct IInterpolater {
  virtual double GetValue(double x[3]) const = 0;
};
}  // namespace mr