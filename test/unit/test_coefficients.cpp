#include "materials.hpp"
#include <gtest/gtest.h>

TEST(CoefficientsTest, CheckData) {
  hephaestus::Subdomain wire("wire", 1);
  wire.property_map.Register("property_one", new mfem::ConstantCoefficient(1.0),
                             true);
  wire.property_map.Register("property_two",
                             new mfem::ConstantCoefficient(150.0), true);
  hephaestus::Subdomain air("air", 2);
  air.property_map.Register("property_one", new mfem::ConstantCoefficient(26.0),
                            true);
  air.property_map.Register("property_two",
                            new mfem::ConstantCoefficient(152.0), true);

  hephaestus::Coefficients domain_properties(
      std::vector<hephaestus::Subdomain>({wire, air}));

  // Verify predefined values
  mfem::IsoparametricTransformation T;
  mfem::IntegrationPoint ip;

  mfem::Coefficient *pw =
      domain_properties.scalar_property_map.Get("property_one");
  T.Attribute = 1;
  EXPECT_FLOAT_EQ(pw->Eval(T, ip), 1.0);
  T.Attribute = 2;
  EXPECT_FLOAT_EQ(pw->Eval(T, ip), 26.0);

  pw = domain_properties.scalar_property_map.Get("property_two");
  T.Attribute = 1;
  EXPECT_FLOAT_EQ(pw->Eval(T, ip), 150.0);
  T.Attribute = 2;
  EXPECT_FLOAT_EQ(pw->Eval(T, ip), 152.0);
}
