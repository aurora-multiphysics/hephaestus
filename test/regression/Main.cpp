#include "mfem.hpp"
#include <iostream>
#include <catch2/catch_session.hpp>

const char *DATA_DIR = "../data/";

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  int result = Catch::Session().run( argc, argv );
  MPI_Finalize();
  
  return result;
}
