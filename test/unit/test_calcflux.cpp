#include "open_coil.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

static void JExact(const mfem::Vector & x, mfem::Vector & J);

extern const char * DATA_DIR;

// Floating point error tolerance
const double eps{1e-10};

TEST_CASE("Flux calculation", "[CalcFlux]")
{

  int par_ref_lvl = -1;
  int order = 1;

  mfem::Mesh mesh((std::string(DATA_DIR) + "team7.g").c_str(), 1, 1);

  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));

  for (int l = 0; l < par_ref_lvl; ++l)
    pmesh->UniformRefinement();

  mfem::ND_FECollection HCurl_Collection(order, pmesh.get()->Dimension());
  mfem::ParFiniteElementSpace HCurlFESpace(pmesh.get(), &HCurl_Collection);
  mfem::ParGridFunction j(&HCurlFESpace);

  mfem::VectorFunctionCoefficient vecField(3, JExact);
  j.ProjectCoefficient(vecField);

  double area = 0.0025;
  double theta = M_PI / 4.;

  double flux = hephaestus::calcFlux(&j, 7);

  REQUIRE_THAT(flux, Catch::Matchers::WithinAbs(area * cos(theta), eps));
}

static void
JExact(const mfem::Vector & x, mfem::Vector & J)
{
  J = 0.0;
  J(0) = 1.0 / sqrt(2);
  J(1) = 1.0 / sqrt(2);
}
