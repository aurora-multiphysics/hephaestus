#include "mfem.hpp"
#include <gtest/gtest.h>
#include <iostream>

const char *DATA_DIR = "../data/";

int main(int argc, char *argv[]) {
  std::cout << "Unit tests from test/unit/Main.cpp" << std::endl;
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&DATA_DIR, "-dataDir", "--data_directory",
                 "Directory storing input data for tests.");
  args.Parse();
  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
