#include "joule_solver.hpp"
#include <gtest/gtest.h>

TEST(dummyTest, CheckData) {
  EXPECT_EQ(1.0e-9, mfem::electromagnetics::SOLVER_TOL);
}
