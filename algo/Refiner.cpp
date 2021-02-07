#include "Refiner.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPointInterpolator.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>

namespace mr {

static constexpr double EPISILON = 0.001;

Refiner::Refiner(vtkPolyData* high_res_points, RefinerConfig config)
    : high_res_points_{high_res_points}, config_{config} {
  locator_->SetDataSet(high_res_points);
  locator_->BuildLocator();
}

vtkNew<vtkPolyData> Refiner::Refine(vtkPolyData* mesh) const {
  vtkPolyData* input_ds = vtkPolyData::New();

  input_ds->CopyStructure(mesh);
  input_ds->GetPointData()->PassData(mesh->GetPointData());
  input_ds->GetCellData()->PassData(mesh->GetCellData());

  for (int level = 0; level < config_.refine_times; ++level) {
    input_ds->BuildLinks();
    vtkIdType n_cells = input_ds->GetNumberOfCells();

    vtkNew<vtkIdList> point_ids_to_refine = FindPointsToRefine(input_ds);
    vtkNew<vtkIdList> cell_ids_to_refine =
        FindUniqueCellsByPoints(input_ds, point_ids_to_refine);

    vtkNew<vtkPoints> output_points;
    output_points->SetDataTypeToDouble();
    output_points->GetData()->DeepCopy(input_ds->GetPoints()->GetData());

    vtkNew<vtkPointData> output_pd;
    output_pd->CopyAllocate(input_ds->GetPointData(),
                            2 * input_ds->GetNumberOfPoints());

    vtkNew<vtkCellData> output_cd;
    output_cd->CopyAllocate(input_ds->GetCellData(),
                            n_cells + 3 * cell_ids_to_refine->GetNumberOfIds());

    vtkNew<vtkCellArray> output_cells;
    output_cells->AllocateEstimate(
        n_cells + 3 * cell_ids_to_refine->GetNumberOfIds(), 3);

    vtkNew<vtkEdgeTable> edge_table;
    edge_table->InitEdgeInsertion(input_ds->GetNumberOfPoints());

    vtkNew<vtkIntArray> edge_data;
    edge_data->SetNumberOfComponents(3);
    edge_data->SetNumberOfTuples(n_cells);

    GenerateSubdivisionPoints(input_ds, edge_table, edge_data, output_points,
                              output_pd, cell_ids_to_refine);
    GenerateSubdivisionCells(input_ds, edge_table, edge_data, output_cells,
                             output_cd);

    input_ds->Delete();
    input_ds = vtkPolyData::New();
    input_ds->SetPoints(output_points);
    input_ds->SetPolys(output_cells);
    input_ds->GetPointData()->PassData(output_pd);
    input_ds->GetCellData()->PassData(output_cd);
    input_ds->Squeeze();
  }

  vtkNew<vtkPolyData> res;
  res->SetPoints(input_ds->GetPoints());
  res->SetPolys(input_ds->GetPolys());
  res->GetPointData()->PassData(input_ds->GetPointData());
  res->GetCellData()->PassData(input_ds->GetCellData());
  input_ds->Delete();

  return res;
}

vtkIdType Refiner::FindEdge(vtkPolyData* mesh, vtkIdType cell_id, vtkIdType p1,
                            vtkIdType p2, vtkIntArray* edge_data,
                            vtkIdList* cell_ids) {
  mesh->GetCellEdgeNeighbors(cell_id, p1, p2, cell_ids);

  for (vtkIdType i = 0; i < cell_ids->GetNumberOfIds(); ++i) {
    vtkIdType current_cell_id = cell_ids->GetId(i);
    vtkCell* cell = mesh->GetCell(current_cell_id);
    int n_edges = cell->GetNumberOfEdges();

    for (int edge_id = 0; edge_id < n_edges; ++edge_id) {
      vtkIdType tp1 = cell->GetPointId((2 + edge_id) % 3);
      vtkIdType tp2 = cell->GetPointId(edge_id);

      if ((tp1 == p1 && tp2 == p2) || (tp1 == p2 && tp2 == p1)) {
        return (int)edge_data->GetComponent(current_cell_id, edge_id);
      }
    }
  }

  throw std::runtime_error("Edge should be found");
}

vtkIdType Refiner::InterpolatePosition(vtkPoints* input_points,
                                       vtkPoints* output_points,
                                       vtkIdList* stencil,
                                       double* weights) const {
  double x[3] = {0., 0., 0.};

  for (vtkIdType i = 0; i < stencil->GetNumberOfIds(); ++i) {
    double input_x[3];
    input_points->GetPoint(stencil->GetId(i), input_x);
    for (int j = 0; j < 3; ++j) {
      x[j] += input_x[j] * weights[i];
    }
  }

  double expected_z = GetInterpolatedZFromHighResData(x);
  if (expected_z > EPISILON) {
    assert(abs(expected_z - x[2]) < 20);

    x[2] = (x[2] * config_.original_z_weight +
            expected_z * (1 - config_.original_z_weight));
  }

  return output_points->InsertNextPoint(x);
}

double Refiner::GetInterpolatedZFromHighResData(double* x) const {
  double expected_z = 0;
  double sum = 0.;
  vtkNew<vtkIdList> point_ids_in_radius;
  locator_->FindPointsWithinRadius(config_.sample_radius, x,
                                   point_ids_in_radius);
  if (point_ids_in_radius->GetNumberOfIds() > 0) {
    for (auto point_id_in_radius : *point_ids_in_radius) {
      double p[3];
      high_res_points_->GetPoint(point_id_in_radius, p);
      double d =
          sqrt((p[0] - x[0]) * (p[0] - x[0]) + (p[1] - x[1]) * (p[1] - x[1]));

      if (d < EPISILON) continue;

      double d_p = pow(d, config_.idw_p);

      expected_z += p[2] / d_p;
      sum += 1. / d_p;
    }

    expected_z /= sum;
  }

  return expected_z;
}

void Refiner::GenerateSubdivisionPoints(vtkPolyData* input_ds,
                                        vtkEdgeTable* edge_table,
                                        vtkIntArray* edge_data,
                                        vtkPoints* output_points,
                                        vtkPointData* output_pd,
                                        vtkIdList* cell_ids_to_refine) const {
  static double weights[2] = {.5, .5};
  vtkIdType n_cells = input_ds->GetNumberOfCells();

  // To avoid frequent memory allocation
  vtkNew<vtkIdList> temp_cell_ids;

  vtkNew<vtkIdList> temp_point_ids;
  temp_point_ids->SetNumberOfIds(2);

  for (vtkIdType cell_id : *cell_ids_to_refine) {
    const vtkIdType* cell_point_ids;
    input_ds->GetCell(cell_id, cell_point_ids);

    for (int i = 0; i < 3; ++i) {
      vtkIdType new_id;

      vtkIdType p1 = cell_point_ids[(i + 2) % 3];
      vtkIdType p2 = cell_point_ids[i];

      if (edge_table->IsEdge(p1, p2) == -1) {
        edge_table->InsertEdge(p1, p2);

        temp_point_ids->SetId(0, p1);
        temp_point_ids->SetId(1, p2);

        new_id = InterpolatePosition(input_ds->GetPoints(), output_points,
                                     temp_point_ids, weights);
        output_pd->InterpolatePoint(input_ds->GetPointData(), new_id,
                                    temp_point_ids, weights);
      } else {
        new_id = FindEdge(input_ds, cell_id, p1, p2, edge_data, temp_cell_ids);
      }

      edge_data->InsertComponent(cell_id, i, new_id);
    }
  }

  for (vtkIdType cell_id = 0; cell_id < n_cells; ++cell_id) {
    if (cell_ids_to_refine->IsId(cell_id) != -1) continue;

    const vtkIdType* cell_point_ids;
    input_ds->GetCell(cell_id, cell_point_ids);

    for (int i = 0; i < 3; ++i) {
      vtkIdType new_id;

      vtkIdType p1 = cell_point_ids[(i + 2) % 3];
      vtkIdType p2 = cell_point_ids[i];
      if (edge_table->IsEdge(p1, p2) != -1) {
        new_id = FindEdge(input_ds, cell_id, p1, p2, edge_data, temp_cell_ids);
        edge_data->InsertComponent(cell_id, i, new_id);
      }
    }
  }
}

void Refiner::GenerateSubdivisionCells(vtkPolyData* input_ds,
                                       vtkEdgeTable* edge_table,
                                       vtkIntArray* edge_data,
                                       vtkCellArray* output_cells,
                                       vtkCellData* output_cd) const {
  vtkIdType n_cells = input_ds->GetNumberOfCells();

  vtkCellData* input_cd = input_ds->GetCellData();

  for (vtkIdType cell_id = 0; cell_id < n_cells; ++cell_id) {
    vtkIdType n_points;
    const vtkIdType* point_ids;
    input_ds->GetCellPoints(cell_id, n_points, point_ids);
    if (n_points != 3) throw std::runtime_error("Only support triangles");

    int n_interpolated_edges = 0;
    int edge_points[3] = {-1, -1, -1};
    for (int i = 0; i < 3; ++i) {
      vtkIdType p1 = point_ids[(i + 2) % 3];
      vtkIdType p2 = point_ids[i];

      if (edge_table->IsEdge(p1, p2) != -1) {
        ++n_interpolated_edges;
        edge_points[i] = (int)edge_data->GetComponent(cell_id, i);
      }
    }

    switch (n_interpolated_edges) {
      case 0: {
        vtkIdType new_id = output_cells->InsertNextCell(3, point_ids);
        output_cd->CopyData(input_cd, cell_id, new_id);
      } break;
      case 1: {
        vtkIdType new_points_ids[3];
        vtkIdType new_id;
        for (int i = 0; i < 3; ++i) {
          if (edge_points[i] == -1) continue;

          new_points_ids[0] = point_ids[(i + 2) % 3];
          new_points_ids[1] = point_ids[(i + 1) % 3];
          new_points_ids[2] = edge_points[i];
          new_id = output_cells->InsertNextCell(3, new_points_ids);
          output_cd->CopyData(input_cd, cell_id, new_id);

          new_points_ids[0] = point_ids[i];
          new_points_ids[1] = point_ids[(i + 1) % 3];
          new_points_ids[2] = edge_points[i];
          new_id = output_cells->InsertNextCell(3, new_points_ids);
          output_cd->CopyData(input_cd, cell_id, new_id);
        }
      } break;
      case 2: {
        vtkIdType new_points_ids[3];
        vtkIdType new_id;
        for (int i = 0; i < 3; ++i) {
          if (edge_points[i] != -1) continue;

          vtkIdType p1 = (i + 2) % 3;
          vtkIdType p2 = i;
          vtkIdType p3 = (i + 1) % 3;

          new_points_ids[0] = point_ids[p3];
          new_points_ids[1] = edge_points[(i + 1) % 3];
          new_points_ids[2] = edge_points[(i + 2) % 3];
          new_id = output_cells->InsertNextCell(3, new_points_ids);
          output_cd->CopyData(input_cd, cell_id, new_id);

          new_points_ids[0] = point_ids[p2];
          new_points_ids[1] = edge_points[(i + 2) % 3];
          new_points_ids[2] = edge_points[(i + 1) % 3];
          new_id = output_cells->InsertNextCell(3, new_points_ids);
          output_cd->CopyData(input_cd, cell_id, new_id);

          new_points_ids[0] = point_ids[p2];
          new_points_ids[1] = point_ids[p1];
          new_points_ids[2] = edge_points[(i + 2) % 3];
          new_id = output_cells->InsertNextCell(3, new_points_ids);
          output_cd->CopyData(input_cd, cell_id, new_id);
        }
      } break;
      case 3: {
        vtkIdType new_points_ids[3];
        for (int i = 0; i < 3; ++i) {
          new_points_ids[i] = point_ids[i];
          new_points_ids[(i + 1) % 3] = edge_points[(1 + i) % 3];
          new_points_ids[(i + 2) % 3] = edge_points[i];
          vtkIdType new_id = output_cells->InsertNextCell(3, new_points_ids);
          output_cd->CopyData(input_cd, cell_id, new_id);
        }

        new_points_ids[0] = edge_points[1];
        new_points_ids[1] = edge_points[2];
        new_points_ids[2] = edge_points[0];
        vtkIdType new_id = output_cells->InsertNextCell(3, new_points_ids);
        output_cd->CopyData(input_cd, cell_id, new_id);
      } break;
      default:
        throw std::runtime_error(
            "The number of interpolated edges shouldn't be greater than 3");
    }
  }
}

vtkNew<vtkIdList> Refiner::FindUniqueCellsByPoints(vtkPolyData* mesh,
                                                   vtkIdList* point_ids) const {
  vtkNew<vtkIdList> res;

  for (vtkIdType point_id : *point_ids) {
    vtkNew<vtkIdList> cell_ids;
    mesh->GetPointCells(point_id, cell_ids);

    for (vtkIdType j = 0; j < cell_ids->GetNumberOfIds(); ++j) {
      vtkIdType cell_id = cell_ids->GetId(j);

      vtkIdType n_points;
      const vtkIdType* cell_point_ids;
      mesh->GetCellPoints(cell_id, n_points, cell_point_ids);
      if (n_points != 3) throw std::runtime_error("Only support triangles");

      double p[3][3];
      for (int i = 0; i < 3; ++i) mesh->GetPoint(cell_point_ids[i], p[i]);
      auto area = vtkTriangle::TriangleArea(p[0], p[1], p[2]);

      if (area >= config_.min_triangle_area) {
        res->InsertUniqueId(cell_id);
      }
    }
  }

  return res;
}

vtkNew<vtkIdList> Refiner::FindPointsToRefine(vtkPolyData* mesh) const {
  vtkNew<vtkIdList> res;

  auto points = mesh->GetPoints();

  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
    double p[3];
    points->GetPoint(i, p);

    vtkNew<vtkIdList> point_ids_in_radius;
    locator_->FindPointsWithinRadius(config_.judge_radius, p,
                                     point_ids_in_radius);

    for (vtkIdType j = 0; j < point_ids_in_radius->GetNumberOfIds(); ++j) {
      double p_in_radius[3];
      vtkIdType point_in_radius = point_ids_in_radius->GetId(j);
      high_res_points_->GetPoint(point_in_radius, p_in_radius);

      if (abs(p[2] - p_in_radius[2]) >= config_.delta_z_threshold) {
        res->InsertNextId(i);
        break;
      }
    }
  }

  return res;
}

}  // namespace mr
