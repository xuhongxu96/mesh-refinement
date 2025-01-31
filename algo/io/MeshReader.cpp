#include "MeshReader.h"

#include <vtkCellData.h>

#include <fstream>
#include <sstream>
#include <string>

namespace mr {
vtkNew<vtkPolyData> MeshReader::Read(std::istream& is) const {
  using namespace std::string_literals;

  vtkNew<vtkPolyData> res;

  // read points
  {
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();

    uint32_t point_count;
    std::string utm_mode;
    is >> point_count >> utm_mode;

    for (uint32_t i = 0; i < point_count; ++i) {
      uint32_t id;
      double x, y, z;
      int code;  // 0 for internal nodes, 1 for a node on land boundaries, and
                 // >1 for all other boudanries

      is >> id >> x >> y >> z >> code;

      if (id != i + 1) {
        throw std::runtime_error("Point id in mesh file should be in order: "s +
                                 std::to_string(id));
      }

      points->InsertNextPoint(x, y, z);
    }

    res->SetPoints(points);
  }

  // read strips
  {
    vtkNew<vtkCellArray> strips;
    uint32_t strip_count;
    int max_nodes, strip_type;
    is >> strip_count >> max_nodes >> strip_type;

    if (max_nodes != 3 || strip_type != 21) {
      throw std::runtime_error("Only triangle strips are supported");
    }

    for (uint32_t i = 0; i < strip_count; ++i) {
      uint32_t id;
      uint32_t p0, p1, p2;

      is >> id >> p0 >> p1 >> p2;

      if (id != i + 1) {
        throw std::runtime_error("Strip id in mesh file should be in order: "s +
                                 std::to_string(id));
      }

      strips->InsertNextCell(3);
      strips->InsertCellPoint(p0 - 1);
      strips->InsertCellPoint(p1 - 1);
      strips->InsertCellPoint(p2 - 1);

      res->SetPolys(strips);
    }
  }

  return res;
}

vtkNew<vtkPolyData> MeshReader::Read(const std::string& path) const {
  std::ifstream ifs(path);
  return Read(ifs);
}

}  // namespace mr
