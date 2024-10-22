#include <gtest/gtest.h>
#include <io/MeshReader.h>
TEST(MeshReader, Test) {
  mr::MeshReader reader;
  auto mesh = reader.Read(DATA_PATH "original.mesh");

  ASSERT_EQ(mesh->GetNumberOfPoints(), 18553);
  ASSERT_EQ(mesh->GetNumberOfCells(), 36615);
}