#include "coefficients.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Floating point error tolerance
const double eps{1e-10};

TEST_CASE("CoefficientsTest", "[CheckData]")
{
  hephaestus::Subdomain wire("wire", 1);
  wire._scalar_coefficients.Register("property_one",
                                     std::make_shared<mfem::ConstantCoefficient>(1.0));
  wire._scalar_coefficients.Register("property_two",
                                     std::make_shared<mfem::ConstantCoefficient>(150.0));

  hephaestus::Subdomain air("air", 2);
  air._scalar_coefficients.Register("property_one",
                                    std::make_shared<mfem::ConstantCoefficient>(26.0));
  air._scalar_coefficients.Register("property_two",
                                    std::make_shared<mfem::ConstantCoefficient>(152.0));

  hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({wire, air}));

  // Verify predefined values
  mfem::IsoparametricTransformation t;
  mfem::IntegrationPoint ip;

  mfem::Coefficient * pw = coefficients._scalars.GetPtr("property_one");
  t.Attribute = 1;
  REQUIRE_THAT(pw->Eval(t, ip), Catch::Matchers::WithinAbs(1.0, eps));
  t.Attribute = 2;
  REQUIRE_THAT(pw->Eval(t, ip), Catch::Matchers::WithinAbs(26.0, eps));

  pw = coefficients._scalars.GetPtr("property_two");
  t.Attribute = 1;
  REQUIRE_THAT(pw->Eval(t, ip), Catch::Matchers::WithinAbs(150.0, eps));
  t.Attribute = 2;
  REQUIRE_THAT(pw->Eval(t, ip), Catch::Matchers::WithinAbs(152.0, eps));
}
