#pragma once

#include <vtkPolyData.h>

#include <iostream>
#include <string>

namespace mr {

class MeshWriter {
 public:
  void Write(std::ostream& os, vtkPolyData* mesh) const;
  void Write(const std::string& path, vtkPolyData* mesh) const;
};

}  // namespace mr