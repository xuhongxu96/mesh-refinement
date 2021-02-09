#include "MeshWriter.h"

#include <vtkFeatureEdges.h>
#include <vtkIdFilter.h>
#include <vtkPointData.h>

#include <fstream>
#include <unordered_set>

namespace mr {

static constexpr auto delimiter = " ";

void MeshWriter::Write(std::ostream& os, vtkPolyData* mesh) const {
  auto n_points = mesh->GetNumberOfPoints();
  auto n_cells = mesh->GetNumberOfCells();

  vtkNew<vtkIdFilter> id_filter;
  id_filter->SetInputData(mesh);
  id_filter->SetPointIdsArrayName("ids");
  id_filter->SetCellIds(false);
  id_filter->Update();

  vtkNew<vtkFeatureEdges> feature_edges;
  feature_edges->SetInputConnection(id_filter->GetOutputPort());
  feature_edges->BoundaryEdgesOn();
  feature_edges->FeatureEdgesOff();
  feature_edges->ManifoldEdgesOff();
  feature_edges->NonManifoldEdgesOff();
  feature_edges->Update();

  auto boundary_arr =
      feature_edges->GetOutput()->GetPointData()->GetArray("ids");

  std::unordered_set<vtkIdType> boundary_point_id_set;
  for (vtkIdType i = 0; i < feature_edges->GetOutput()->GetNumberOfPoints();
       ++i) {
    boundary_point_id_set.insert(boundary_arr->GetTuple1(i));
  }

  os.precision(8);
  os.setf(std::ios::fixed);

  os << n_points << delimiter << "NON-UTM" << std::endl;

  for (vtkIdType i = 0; i < n_points; ++i) {
    double p[3];
    mesh->GetPoint(i, p);

    int boundary_code = 0;
    if (boundary_point_id_set.find(i) != boundary_point_id_set.end()) {
      boundary_code = 1;
    }

    os << i + 1 << delimiter << p[0] << delimiter << p[1] << delimiter << p[2]
       << delimiter << boundary_code << std::endl;
  }

  os << n_cells << delimiter << 3 << delimiter << 21 << delimiter << std::endl;

  for (vtkIdType i = 0; i < n_cells; ++i) {
    vtkIdType n_pt;
    const vtkIdType* point_ids;
    mesh->GetCellPoints(i, n_pt, point_ids);
    if (n_pt != 3) throw std::runtime_error("Only support triangles");

    os << i + 1 << delimiter << point_ids[0] + 1 << delimiter
       << point_ids[1] + 1 << delimiter << point_ids[2] + 1 << std::endl;
  }
}

void MeshWriter::Write(const std::string& path, vtkPolyData* mesh) const {
  std::ofstream ofs(path);
  Write(ofs, mesh);
}

}  // namespace mr