#include <gtest/gtest.h>
#include "mfem.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << "Integration tests from test/integration/Main.cpp" << std::endl;
    int result = 0;

    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    result = RUN_ALL_TESTS();
    MPI_Finalize();

    return result;

}
