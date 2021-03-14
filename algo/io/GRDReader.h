#pragma once

#include <vtkPolyData.h>

#include <iostream>
#include <string>

namespace mr {

class GRDReader {
 public:
  vtkNew<vtkPolyData> Read(std::istream& is) const;
  vtkNew<vtkPolyData> Read(const std::string& path) const;
};

}  // namespace mr