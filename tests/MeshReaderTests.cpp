#include <MeshReader.h>
#include <gtest/gtest.h>
TEST(MeshReader, Test) {
  mr::MeshReader reader;
  mr::Mesh mesh = reader.Read(DATA_PATH "original.mesh");

  ASSERT_EQ(mesh.points->GetNumberOfPoints(), 18553);
  ASSERT_EQ(mesh.strips->GetNumberOfCells(), 36615);
  ASSERT_EQ(mesh.poly_data->GetPoints(), mesh.points);
  ASSERT_EQ(mesh.poly_data->GetStrips(), mesh.strips);
}