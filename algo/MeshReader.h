#pragma once

#include <vtkPoints.h>
#include <vtkPolyData.h>

#include <string_view>

namespace mr {

struct Mesh {
  vtkNew<vtkPoints> points;
  vtkNew<vtkCellArray> strips;
  vtkNew<vtkPolyData> poly_data;
};

class MeshReader {
 public:
  Mesh Read(const std::string& path) const;
};

}  // namespace mr