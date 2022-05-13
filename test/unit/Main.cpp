#include "mfem.hpp"
#include <gtest/gtest.h>
#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << "Unit tests from test/unit/Main.cpp" << std::endl;
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
