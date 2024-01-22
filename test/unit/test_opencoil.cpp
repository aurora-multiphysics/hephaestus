#include "open_coil.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

extern const char * DATA_DIR;

TEST_CASE("OpenCoilTest", "[CheckData]")
{

  // Floating point error tolerance
  const double eps{1e-10};

  int par_ref_lvl = -1;
  int order = 1;

  mfem::Mesh mesh((std::string(DATA_DIR) + "coil.gen").c_str(), 1, 1);
  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

  for (int l = 0; l < par_ref_lvl; ++l)
    pmesh->UniformRefinement();

  mfem::ND_FECollection h_curl_collection(order, pmesh.get()->Dimension());
  mfem::ParFiniteElementSpace h_curl_fe_space(pmesh.get(), &h_curl_collection);
  mfem::ParGridFunction e(&h_curl_fe_space);

  const double ival = 10.0;
  const double cond_val = 1e6;

  mfem::ConstantCoefficient itot(ival);
  mfem::ConstantCoefficient conductivity(cond_val);

  hephaestus::InputParameters ocs_params;
  hephaestus::BCMap bc_maps;

  hephaestus::Coefficients coefficients;
  coefficients._scalars.Register(std::string("Itotal"), &itot, false);
  coefficients._scalars.Register(std::string("Conductivity"), &conductivity, false);

  hephaestus::FESpaces fespaces;
  fespaces.Register(std::string("HCurl"), &h_curl_fe_space, false);

  hephaestus::GridFunctions gridfunctions;
  gridfunctions.Register(std::string("E"), &e, false);

  ocs_params.SetParam("GradPotentialName", std::string("E"));
  ocs_params.SetParam("IFuncCoefName", std::string("Itotal"));
  ocs_params.SetParam("PotentialName", std::string("V"));
  ocs_params.SetParam("ConductivityCoefName", std::string("Conductivity"));

  std::pair<int, int> elec_attrs{1, 2};
  mfem::Array<int> submesh_domains;
  submesh_domains.Append(1);

  hephaestus::OpenCoilSolver opencoil(ocs_params, submesh_domains, elec_attrs);
  opencoil.Init(gridfunctions, fespaces, bc_maps, coefficients);
  mfem::ParLinearForm dummy(&h_curl_fe_space);
  opencoil.Apply(&dummy);

  double flux = hephaestus::calcFlux(&e, elec_attrs.first, conductivity);

  REQUIRE_THAT(flux, Catch::Matchers::WithinAbs(ival, eps));
}
