#pragma once

#include <vtkPolyDataAlgorithm.h>

class vtkAbstractTransform;
class vtkCellArray;
class vtkIdList;
class vtkPointSet;

#define VTK_DELAUNAY_XY_PLANE 0
#define VTK_SET_TRANSFORM_PLANE 1
#define VTK_BEST_FITTING_PLANE 2

class xDelaunay2D : public vtkPolyDataAlgorithm {
 public:
  vtkTypeMacro(xDelaunay2D, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Construct object with Alpha = 0.0; Tolerance = 0.001; Offset = 1.25;
   * BoundingTriangulation turned off.
   */
  static xDelaunay2D* New();

  /**
   * Specify the source object used to specify constrained edges and loops.
   * (This is optional.) If set, and lines/polygons are defined, a constrained
   * triangulation is created. The lines/polygons are assumed to reference
   * points in the input point set (i.e. point ids are identical in the
   * input and source).
   * Note that this method does not connect the pipeline. See
   * SetSourceConnection for connecting the pipeline.
   */
  void SetSourceData(vtkPolyData*);

  /**
   * Specify the source object used to specify constrained edges and loops.
   * (This is optional.) If set, and lines/polygons are defined, a constrained
   * triangulation is created. The lines/polygons are assumed to reference
   * points in the input point set (i.e. point ids are identical in the
   * input and source).
   * New style. This method is equivalent to SetInputConnection(1, algOutput).
   */
  void SetSourceConnection(vtkAlgorithmOutput* algOutput);

  /**
   * Get a pointer to the source object.
   */
  vtkPolyData* GetSource();

  //@{
  /**
   * Specify alpha (or distance) value to control output of this filter.
   * For a non-zero alpha value, only edges or triangles contained within
   * a sphere centered at mesh vertices will be output. Otherwise, only
   * triangles will be output.
   */
  vtkSetClampMacro(Alpha, double, 0.0, VTK_DOUBLE_MAX);
  vtkGetMacro(Alpha, double);
  //@}

  //@{
  /**
   * Specify a tolerance to control discarding of closely spaced points.
   * This tolerance is specified as a fraction of the diagonal length of
   * the bounding box of the points.
   */
  vtkSetClampMacro(Tolerance, double, 0.0, 1.0);
  vtkGetMacro(Tolerance, double);
  //@}

  //@{
  /**
   * Specify a multiplier to control the size of the initial, bounding
   * Delaunay triangulation.
   */
  vtkSetClampMacro(Offset, double, 0.75, VTK_DOUBLE_MAX);
  vtkGetMacro(Offset, double);
  //@}

  //@{
  /**
   * Boolean controls whether bounding triangulation points (and associated
   * triangles) are included in the output. (These are introduced as an
   * initial triangulation to begin the triangulation process. This feature
   * is nice for debugging output.)
   */
  vtkSetMacro(BoundingTriangulation, vtkTypeBool);
  vtkGetMacro(BoundingTriangulation, vtkTypeBool);
  vtkBooleanMacro(BoundingTriangulation, vtkTypeBool);
  //@}

  //@{
  /**
   * Set / get the transform which is applied to points to generate a
   * 2D problem.  This maps a 3D dataset into a 2D dataset where
   * triangulation can be done on the XY plane.  The points are
   * transformed and triangulated.  The topology of triangulated
   * points is used as the output topology.  The output points are the
   * original (untransformed) points.  The transform can be any
   * subclass of vtkAbstractTransform (thus it does not need to be a
   * linear or invertible transform).
   */
  virtual void SetTransform(vtkAbstractTransform*);
  vtkGetObjectMacro(Transform, vtkAbstractTransform);
  //@}

  //@{
  /**
   * Define the method to project the input 3D points into a 2D plane for
   * triangulation. When the VTK_DELAUNAY_XY_PLANE is set, the z-coordinate
   * is simply ignored. When VTK_SET_TRANSFORM_PLANE is set, then a transform
   * must be supplied and the points are transformed using it. Finally, if
   * VTK_BEST_FITTING_PLANE is set, then the filter computes a best fitting
   * plane and projects the points onto it.
   */
  vtkSetClampMacro(ProjectionPlaneMode, int, VTK_DELAUNAY_XY_PLANE,
                   VTK_BEST_FITTING_PLANE);
  vtkGetMacro(ProjectionPlaneMode, int);
  //@}

  /**
   * This method computes the best fit plane to a set of points represented
   * by a vtkPointSet. The method constructs a transform and returns it on
   * successful completion (null otherwise). The user is responsible for
   * deleting the transform instance.
   */
  static vtkAbstractTransform* ComputeBestFittingPlane(vtkPointSet* input);

 protected:
  xDelaunay2D();
  ~xDelaunay2D() override;

  int RequestData(vtkInformation*, vtkInformationVector**,
                  vtkInformationVector*) override;

  double Alpha;
  double Tolerance;
  vtkTypeBool BoundingTriangulation;
  double Offset;

  vtkAbstractTransform* Transform;

  int ProjectionPlaneMode;  // selects the plane in 3D where the Delaunay
                            // triangulation will be computed.

 private:
  vtkPolyData* Mesh;  // the created mesh
  double* Points;     // the raw points in double precision
  void SetPoint(vtkIdType id, double* x) {
    vtkIdType idx = 3 * id;
    this->Points[idx] = x[0];
    this->Points[idx + 1] = x[1];
    this->Points[idx + 2] = x[2];
  }

  void GetPoint(vtkIdType id, double x[3]) {
    double* ptr = this->Points + 3 * id;
    x[0] = *ptr++;
    x[1] = *ptr++;
    x[2] = *ptr;
  }

  int NumberOfDuplicatePoints;
  int NumberOfDegeneracies;

  int* RecoverBoundary(vtkPolyData* source);
  int RecoverEdge(vtkPolyData* source, vtkIdType p1, vtkIdType p2);
  void FillPolygons(vtkCellArray* polys, int* triUse);

  int InCircle(double x[3], double x1[3], double x2[3], double x3[3]);
  vtkIdType FindTriangle(double x[3], vtkIdType ptIds[3], vtkIdType tri,
                         double tol, vtkIdType nei[3], vtkIdList* neighbors);
  void CheckEdge(vtkIdType ptId, double x[3], vtkIdType p1, vtkIdType p2,
                 vtkIdType tri, bool recursive);

  int FillInputPortInformation(int, vtkInformation*) override;

 private:
  xDelaunay2D(const xDelaunay2D&) = delete;
  void operator=(const xDelaunay2D&) = delete;
};
