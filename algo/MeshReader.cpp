#include "MeshReader.h"

#include <fstream>
#include <sstream>
#include <string>

namespace mr {

Mesh MeshReader::Read(const std::string& path) const {
  using namespace std::string_literals;

  Mesh res;

  std::ifstream ifs(path);

  // read points
  {
    uint32_t point_count;
    std::string utm_mode;
    ifs >> point_count >> utm_mode;

    for (uint32_t i = 0; i < point_count; ++i) {
      uint32_t id;
      double x, y, z;
      int code;  // 0 for internal nodes, 1 for a node on land boundaries, and
                 // >1 for all other boudanries

      ifs >> id >> x >> y >> z >> code;

      if (id != i + 1) {
        throw std::runtime_error("Point id in mesh file should be in order: "s +
                                 std::to_string(id));
      }

      res.points->InsertNextPoint(x, y, z);
    }
  }

  // read strips
  {
    uint32_t strip_count;
    int max_nodes, strip_type;
    ifs >> strip_count >> max_nodes >> strip_type;

    if (max_nodes != 3 || strip_type != 21) {
      throw std::runtime_error("Only triangle strips are supported");
    }

    for (uint32_t i = 0; i < strip_count; ++i) {
      uint32_t id;
      uint32_t p0, p1, p2;

      ifs >> id >> p0 >> p1 >> p2;

      if (id != i + 1) {
        throw std::runtime_error("Strip id in mesh file should be in order: "s +
                                 std::to_string(id));
      }

      res.strips->InsertNextCell(3);
      res.strips->InsertCellPoint(p0 - 1);
      res.strips->InsertCellPoint(p1 - 1);
      res.strips->InsertCellPoint(p2 - 1);
    }
  }

  res.poly_data->SetPoints(res.points);
  res.poly_data->SetStrips(res.strips);

  return res;
}

}  // namespace mr
