#include "boundary_condition_base.hpp"
#include "boundary_conditions.hpp"
#include <gtest/gtest.h>

TEST(BoundaryConditionTest, CheckData) {
  hephaestus::BCMap bc_map;
  mfem::Array<int> bdr_attrs({1, 2, 3});
  bc_map.Register(
      "tangential_dEdt",
      new hephaestus::BoundaryCondition(std::string("boundary_1"), bdr_attrs),
      true);

  mfem::Array<int> ess_bdr = bc_map.Get("tangential_dEdt")->bdr_attributes;

  for (int i = 0; i < bdr_attrs.Size(); ++i) {
    EXPECT_EQ(bdr_attrs[i], ess_bdr[i])
        << "Arrays x and y differ at index " << i;
  }
}
