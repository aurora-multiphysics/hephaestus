#include "problem_builder.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

TEST_CASE("VariablesTest", "[CheckSetup]")
{
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);
  hephaestus::TimeDomainProblemBuilder * problem_builder =
      new hephaestus::TimeDomainProblemBuilder();

  std::shared_ptr<mfem::ParMesh> pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);
  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P2"));
  problem_builder->AddGridFunction(std::string("vector_potential"), std::string("HCurl"));
  std::unique_ptr<hephaestus::TimeDomainProblem> problem = problem_builder->ReturnProblem();

  mfem::ParGridFunction * stored_gf = problem->gridfunctions.Get("vector_potential");
  mfem::ParFiniteElementSpace * stored_fespace = problem->fespaces.Get("HCurl");

  REQUIRE(stored_fespace->GetVSize() == stored_gf->ParFESpace()->GetVSize());
  REQUIRE(stored_fespace->GetVSize() > 0);
}
