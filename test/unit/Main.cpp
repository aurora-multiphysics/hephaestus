#include <gtest/gtest.h>
#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << "Unit tests from test/unit/Main.cpp" << std::endl;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
