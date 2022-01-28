#include <gtest/gtest.h>
#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << "Regression tests from test/regression/Main.cpp" << std::endl;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
