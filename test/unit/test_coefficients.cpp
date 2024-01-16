#include "coefficients.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Floating point error tolerance
const double eps{1e-10};

TEST_CASE("CoefficientsTest", "[CheckData]")
{
  hephaestus::Subdomain wire("wire", 1);
  wire.scalar_coefficients.Register("property_one", new mfem::ConstantCoefficient(1.0), true);
  wire.scalar_coefficients.Register("property_two", new mfem::ConstantCoefficient(150.0), true);
  hephaestus::Subdomain air("air", 2);
  air.scalar_coefficients.Register("property_one", new mfem::ConstantCoefficient(26.0), true);
  air.scalar_coefficients.Register("property_two", new mfem::ConstantCoefficient(152.0), true);

  hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({wire, air}));

  // Verify predefined values
  mfem::IsoparametricTransformation T;
  mfem::IntegrationPoint ip;

  mfem::Coefficient * pw = coefficients.scalars.Get("property_one");
  T.Attribute = 1;
  REQUIRE_THAT(pw->Eval(T, ip), Catch::Matchers::WithinAbs(1.0, eps));
  T.Attribute = 2;
  REQUIRE_THAT(pw->Eval(T, ip), Catch::Matchers::WithinAbs(26.0, eps));

  pw = coefficients.scalars.Get("property_two");
  T.Attribute = 1;
  REQUIRE_THAT(pw->Eval(T, ip), Catch::Matchers::WithinAbs(150.0, eps));
  T.Attribute = 2;
  REQUIRE_THAT(pw->Eval(T, ip), Catch::Matchers::WithinAbs(152.0, eps));
}
