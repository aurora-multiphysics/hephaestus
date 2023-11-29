#include "mfem.hpp"
#include <iostream>
#include <catch2/catch_session.hpp>

const char *DATA_DIR = "../data/";

int main(int argc, char *argv[]) {
  std::cout << "Regression tests from test/regression/Main.cpp" << std::endl;
  int result = 0;
  MPI_Init(&argc, &argv);
  result = Catch::Session().run( argc, argv );
  MPI_Finalize();
  return result;
}
