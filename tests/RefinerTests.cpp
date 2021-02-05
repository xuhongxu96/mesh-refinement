#include <MeshReader.h>
#include <Refiner.h>
#include <gtest/gtest.h>
#include <vtkPolyDataMapper.h>
#include <vtkSimplePointsReader.h>

TEST(Refiner, Test) {
  vtkNew<vtkSimplePointsReader> reader;
  reader->SetFileName(DATA_PATH "vectors.txt");
  reader->Update();

  mr::MeshReader mesh_reader;
  mr::Mesh mesh = mesh_reader.Read(DATA_PATH "original.mesh");

  mr::Refiner refiner(reader->GetOutput());
  refiner.Refine(mesh.poly_data);
}