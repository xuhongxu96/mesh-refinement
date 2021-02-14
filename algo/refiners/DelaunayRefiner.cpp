#include "DelaunayRefiner.h"

#include <vtkCellArray.h>
#include <vtkFeatureEdges.h>
#include <vtkIdFilter.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkPolygon.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>

#include <array>
#include <unordered_set>
namespace mr {

constexpr auto VTK_DEL2D_TOLERANCE = 1.0e-014;

// Determine whether point x is inside of circumcircle of triangle
// defined by points (x1, x2, x3). Returns non-zero if inside circle.
// (Note that z-component is ignored.)
static int InCircle(double x[3], double x1[3], double x2[3], double x3[3]) {
  double radius2, center[2], dist2;

  radius2 = vtkTriangle::Circumcircle(x1, x2, x3, center);

  // check if inside/outside circumcircle
  dist2 = (x[0] - center[0]) * (x[0] - center[0]) +
          (x[1] - center[1]) * (x[1] - center[1]);

  if (dist2 < (0.999999999999 * radius2)) {
    return 1;
  } else {
    return 0;
  }
}

// Recursive method to locate triangle containing point. Starts with arbitrary
// triangle (tri) and "walks" towards it. Influenced by some of Guibas and
// Stolfi's work. Returns id of enclosing triangle, or -1 if no triangle
// found. Also, the array nei[3] is used to communicate info about points
// that lie on triangle edges: nei[0] is neighboring triangle id, and nei[1]
// and nei[2] are the vertices defining the edge.
static vtkIdType FindTriangle(vtkPolyData* mesh, double x[3],
                              vtkIdType ptIds[3], vtkIdType tri, double tol,
                              vtkIdType nei[3], vtkIdList* neighbors,
                              int& n_duplicate_points, int& n_degeneracies) {
  int i, j, ir, ic, inside, i2, i3;
  const vtkIdType* pts;
  vtkIdType npts;
  vtkIdType newNei;
  double p[3][3], n[2], vp[2], vx[2], dp, minProj;

  // get local triangle info
  mesh->GetCellPoints(tri, npts, pts);
  for (i = 0; i < 3; i++) {
    ptIds[i] = pts[i];
    mesh->GetPoints()->GetPoint(ptIds[i], p[i]);
  }

  // Randomization (of find edge neighbora) avoids walking in
  // circles in certain weird cases
  srand(tri);
  ir = rand() % 3;
  // evaluate in/out of each edge
  for (inside = 1, minProj = VTK_DEL2D_TOLERANCE, ic = 0; ic < 3; ic++) {
    i = (ir + ic) % 3;
    i2 = (i + 1) % 3;
    i3 = (i + 2) % 3;

    // create a 2D edge normal to define a "half-space"; evaluate points (i.e.,
    // candidate point and other triangle vertex not on this edge).
    n[0] = -(p[i2][1] - p[i][1]);
    n[1] = p[i2][0] - p[i][0];
    vtkMath::Normalize2D(n);

    // compute local vectors
    for (j = 0; j < 2; j++) {
      vp[j] = p[i3][j] - p[i][j];
      vx[j] = x[j] - p[i][j];
    }

    vtkMath::Normalize2D(vp);
    // check for duplicate point
    if (vtkMath::Normalize2D(vx) <= tol * 0.00001) {
      ++n_duplicate_points;
      return -1;
    }
    // see if two points are in opposite half spaces
    dp = vtkMath::Dot2D(n, vx) * (vtkMath::Dot2D(n, vp) < 0 ? -1.0 : 1.0);
    if (dp < VTK_DEL2D_TOLERANCE) {
      if (dp < minProj)  // track edge most orthogonal to point direction
      {
        inside = 0;
        nei[1] = ptIds[i];
        nei[2] = ptIds[i2];
        minProj = dp;
      }
    }  // outside this edge
  }    // for each edge

  if (inside)  // all edges have tested positive
  {
    nei[0] = (-1);
    return tri;
  }

  else if (!inside && (fabs(minProj) < VTK_DEL2D_TOLERANCE))  // on edge
  {
    mesh->GetCellEdgeNeighbors(tri, nei[1], nei[2], neighbors);
    nei[0] = neighbors->GetId(0);
    return tri;
  }

  else  // walk towards point
  {
    mesh->GetCellEdgeNeighbors(tri, nei[1], nei[2], neighbors);
    if ((newNei = neighbors->GetId(0)) == nei[0]) {
      ++n_degeneracies;
      return -1;
    } else {
      nei[0] = tri;
      return FindTriangle(mesh, x, ptIds, newNei, tol, nei, neighbors,
                          n_duplicate_points, n_degeneracies);
    }
  }
}

// Recursive method checks whether edge is Delaunay, and if not, swaps edge.
// Continues until all edges are Delaunay. Points p1 and p2 form the edge in
// question; x is the coordinates of the inserted point; tri is the current
// triangle id.
static void CheckEdge(vtkPolyData* mesh, vtkIdType ptId, double x[3],
                      vtkIdType p1, vtkIdType p2, vtkIdType tri,
                      bool recursive) {
  int i;
  const vtkIdType* pts;
  vtkIdType npts;
  vtkIdType numNei, nei, p3;
  double x1[3], x2[3], x3[3];
  vtkIdList* neighbors;
  vtkIdType swapTri[3];

  mesh->GetPoints()->GetPoint(p1, x1);
  mesh->GetPoints()->GetPoint(p2, x2);

  neighbors = vtkIdList::New();
  neighbors->Allocate(2);

  mesh->GetCellEdgeNeighbors(tri, p1, p2, neighbors);
  numNei = neighbors->GetNumberOfIds();

  if (numNei > 0)  // i.e., not a boundary edge
  {
    // get neighbor info including opposite point
    nei = neighbors->GetId(0);
    mesh->GetCellPoints(nei, npts, pts);
    for (i = 0; i < 2; i++) {
      if (pts[i] != p1 && pts[i] != p2) {
        break;
      }
    }
    p3 = pts[i];
    mesh->GetPoints()->GetPoint(p3, x3);

    // see whether point is in circumcircle
    if (InCircle(x3, x, x1, x2)) {  // swap diagonal
      mesh->RemoveReferenceToCell(p1, tri);
      mesh->RemoveReferenceToCell(p2, nei);
      mesh->ResizeCellList(ptId, 1);
      mesh->AddReferenceToCell(ptId, nei);
      mesh->ResizeCellList(p3, 1);
      mesh->AddReferenceToCell(p3, tri);

      swapTri[0] = ptId;
      swapTri[1] = p3;
      swapTri[2] = p2;
      mesh->ReplaceCell(tri, 3, swapTri);

      swapTri[0] = ptId;
      swapTri[1] = p1;
      swapTri[2] = p3;
      mesh->ReplaceCell(nei, 3, swapTri);

      if (recursive) {
        // two new edges become suspect
        CheckEdge(mesh, ptId, x, p3, p2, tri, true);
        CheckEdge(mesh, ptId, x, p1, p3, nei, true);
      }
    }  // in circle
  }    // interior edge

  neighbors->Delete();
}

DelaunayRefiner::DelaunayRefiner(std::shared_ptr<IInterpolater> interpolater,
                                 DelaunayRefinerConfig config)
    : interpolater_{interpolater}, config_{std::move(config)} {}

vtkNew<vtkPolyData> DelaunayRefiner::Refine(vtkPolyData* input,
                                            vtkIdList* cell_ids_to_refine,
                                            vtkPoints* degen_points) const {
  vtkNew<vtkPolyData> mesh;

  /*
  vtkNew<vtkFeatureEdges> feature_edges;
  feature_edges->SetInputData(input);
  feature_edges->BoundaryEdgesOn();
  feature_edges->FeatureEdgesOff();
  feature_edges->ManifoldEdgesOff();
  feature_edges->NonManifoldEdgesOff();
  feature_edges->Update();
*/
  mesh->DeepCopy(input);

  vtkIdType n_old_points = mesh->GetNumberOfPoints();

  auto point_cell_map =
      GeneratePoints(input, mesh->GetPoints(), cell_ids_to_refine);
  mesh->BuildLinks();

  vtkIdType n_points = mesh->GetNumberOfPoints();
  auto tol = mesh->GetLength();

  vtkNew<vtkIdList> neighbors;
  neighbors->Allocate(2);

  int n_duplicate = 0, n_degeneracies = 0;
  vtkIdType nei[3];
  vtkIdType tri[4];
  tri[0] = 0;

  vtkNew<vtkPoints> new_points;
  std::unordered_map<vtkIdType, vtkIdType> new_point_map;

  for (vtkIdType point_id = 0; point_id < n_old_points; point_id++) {
    double x[3];
    mesh->GetPoint(point_id, x);
    new_points->InsertNextPoint(x);
  }

  for (vtkIdType point_id = n_old_points; point_id < n_points; point_id++) {
    double x[3];
    vtkIdType pts[3];
    vtkIdType nodes[4][3];

    mesh->GetPoint(point_id, x);

    nei[0] = (-1);  // where we are coming from...nowhere initially
    tri[0] = point_cell_map[point_id];

    if ((tri[0] = FindTriangle(mesh, x, pts, tri[0], tol, nei, neighbors,
                               n_duplicate, n_degeneracies)) >= 0) {
      new_point_map[point_id] = new_points->InsertNextPoint(x);

      if (nei[0] < 0)  // in triangle
      {
        // delete this triangle; create three new triangles
        // first triangle is replaced with one of the new ones
        nodes[0][0] = point_id;
        nodes[0][1] = pts[0];
        nodes[0][2] = pts[1];
        mesh->RemoveReferenceToCell(pts[2], tri[0]);
        mesh->ReplaceCell(tri[0], 3, nodes[0]);
        mesh->ResizeCellList(point_id, 1);
        mesh->AddReferenceToCell(point_id, tri[0]);

        // create two new triangles
        nodes[1][0] = point_id;
        nodes[1][1] = pts[1];
        nodes[1][2] = pts[2];
        tri[1] = mesh->InsertNextLinkedCell(VTK_TRIANGLE, 3, nodes[1]);

        nodes[2][0] = point_id;
        nodes[2][1] = pts[2];
        nodes[2][2] = pts[0];
        tri[2] = mesh->InsertNextLinkedCell(VTK_TRIANGLE, 3, nodes[2]);

        // Check edge neighbors for Delaunay criterion. If not satisfied, flip
        // edge diagonal. (This is done recursively.)
        CheckEdge(mesh, point_id, x, pts[0], pts[1], tri[0], true);
        CheckEdge(mesh, point_id, x, pts[1], pts[2], tri[1], true);
        CheckEdge(mesh, point_id, x, pts[2], pts[0], tri[2], true);
      }

      else  // on triangle edge
      {
        // update cell list
        vtkIdType numNeiPts;
        const vtkIdType* neiPts;
        vtkIdType p1 = 0;
        vtkIdType p2 = 0;
        vtkIdType p3 = 0;

        mesh->GetCellPoints(nei[0], numNeiPts, neiPts);
        for (int i = 0; i < 3; i++) {
          if (neiPts[i] != nei[1] && neiPts[i] != nei[2]) {
            p1 = neiPts[i];
          }
          if (pts[i] != nei[1] && pts[i] != nei[2]) {
            p2 = pts[i];
          }
        }
        mesh->ResizeCellList(p1, 1);
        mesh->ResizeCellList(p2, 1);

        // replace two triangles
        mesh->RemoveReferenceToCell(nei[2], tri[0]);
        mesh->RemoveReferenceToCell(nei[2], nei[0]);
        nodes[0][0] = point_id;
        nodes[0][1] = p2;
        nodes[0][2] = nei[1];
        mesh->ReplaceCell(tri[0], 3, nodes[0]);
        nodes[1][0] = point_id;
        nodes[1][1] = p1;
        nodes[1][2] = nei[1];
        mesh->ReplaceCell(nei[0], 3, nodes[1]);
        mesh->ResizeCellList(point_id, 2);
        mesh->AddReferenceToCell(point_id, tri[0]);
        mesh->AddReferenceToCell(point_id, nei[0]);

        tri[1] = nei[0];

        // create two new triangles
        nodes[2][0] = point_id;
        nodes[2][1] = p2;
        nodes[2][2] = nei[2];
        tri[2] = mesh->InsertNextLinkedCell(VTK_TRIANGLE, 3, nodes[2]);

        nodes[3][0] = point_id;
        nodes[3][1] = p1;
        nodes[3][2] = nei[2];
        tri[3] = mesh->InsertNextLinkedCell(VTK_TRIANGLE, 3, nodes[3]);

        // Check edge neighbors for Delaunay criterion.
        for (int i = 0; i < 4; i++) {
          CheckEdge(mesh, point_id, x, nodes[i][1], nodes[i][2], tri[i], true);
        }
      }
    }  // if triangle found

    else {
      tri[0] = 0;  // no triangle found
      degen_points->InsertNextPoint(x);
    }
  }

  vtkNew<vtkPolyData> res;
  res->SetPoints(new_points);

  vtkNew<vtkCellArray> cells;

  for (vtkIdType cell_id = 0; cell_id < mesh->GetNumberOfCells(); ++cell_id) {
    vtkNew<vtkIdList> ids;
    mesh->GetCellPoints(cell_id, ids);
    for (vtkIdType i = 0; i < ids->GetNumberOfIds(); i++) {
      if (auto it = new_point_map.find(ids->GetId(i));
          it != new_point_map.end()) {
        ids->SetId(i, it->second);
      }
    }
    cells->InsertNextCell(ids);
  }

  res->SetPolys(cells);
  res->BuildLinks();

  return res;
}

std::unordered_map<vtkIdType, vtkIdType> DelaunayRefiner::GeneratePoints(
    vtkPolyData* input, vtkPoints* output,
    vtkIdList* cell_ids_to_refine) const {
  std::unordered_map<vtkIdType, vtkIdType> res;
  for (vtkIdType cell_id : *cell_ids_to_refine) {
    vtkIdType n_points;
    const vtkIdType* point_ids;
    input->GetCellPoints(cell_id, n_points, point_ids);
    assert(n_points == 3);

    double p[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
      double cell_p[3];
      input->GetPoint(point_ids[i], cell_p);
      for (int j = 0; j < 3; ++j) p[j] += cell_p[j];
    }

    for (int j = 0; j < 3; ++j) p[j] /= 3;

    double z = interpolater_->GetValue(p);
    if (z == 0) continue;

    res[output->InsertNextPoint(p[0], p[1], z)] = cell_id;
  }

  return res;
}

}  // namespace mr