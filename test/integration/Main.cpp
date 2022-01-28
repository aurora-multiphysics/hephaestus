#include <gtest/gtest.h>
#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << "Integration tests from test/integration/Main.cpp" << std::endl;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
