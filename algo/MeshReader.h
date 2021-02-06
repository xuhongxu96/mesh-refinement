#pragma once

#include <vtkPolyData.h>

#include <string_view>

namespace mr {

class MeshReader {
 public:
  vtkNew<vtkPolyData> Read(const std::string& path) const;
};

}  // namespace mr