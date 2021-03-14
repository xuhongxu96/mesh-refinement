#include "GRDReader.h"

#include <vtkCellData.h>

#include <fstream>
#include <sstream>
#include <string>

namespace mr {
vtkNew<vtkPolyData> GRDReader::Read(std::istream& is) const {
  using namespace std::string_literals;

  vtkNew<vtkPolyData> res;

  std::string grd_id;
  std::getline(is, grd_id);

  vtkIdType n_cells, n_points;
  is >> n_cells >> n_points;

  // read points
  {
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();

    for (vtkIdType i = 0; i < n_points; ++i) {
      uint32_t id;
      double x, y, z;
      int code;  // 0 for internal nodes, 1 for a node on land boundaries, and
                 // >1 for all other boudanries

      is >> id >> x >> y >> z;

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

    for (vtkIdType i = 0; i < n_cells; ++i) {
      uint32_t id;
      uint32_t n_pts;

      is >> id >> n_pts;

      if (n_pts != 3) throw std::runtime_error("Only support triangles");

      uint32_t p0, p1, p2;
      is >> p0 >> p1 >> p2;

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

vtkNew<vtkPolyData> GRDReader::Read(const std::string& path) const {
  std::ifstream ifs(path);
  return Read(ifs);
}

}  // namespace mr
