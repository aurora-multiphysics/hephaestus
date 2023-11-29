#include "mfem.hpp"
#include <iostream>
#include <catch2/catch_session.hpp>

const char *DATA_DIR = "../data/";

int main(int argc, char *argv[]) {
  std::cout << "Integration tests from test/integration/Main.cpp" << std::endl;
  int result = 0;
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&DATA_DIR, "-dataDir", "--data_directory",
                 "Directory storing input data for tests.");
  args.Parse();
  MPI_Init(&argc, &argv);
  result = Catch::Session().run( argc, argv );
  MPI_Finalize();
  return result;
}
