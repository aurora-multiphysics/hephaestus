#include "problem_builder.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

TEST(VariablesTest, CheckSetup) {
  mfem::Mesh mesh(
      (std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);
  hephaestus::TransientProblemBuilder *problem_builder =
      new hephaestus::TransientProblemBuilder();

  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);
  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P2"));
  problem_builder->AddGridFunction(std::string("vector_potential"),
                                   std::string("HCurl"));
  std::unique_ptr<hephaestus::TransientProblem> problem =
      problem_builder->ReturnProblem();

  mfem::ParGridFunction *stored_gf =
      problem->gridfunctions.Get("vector_potential");
  mfem::ParFiniteElementSpace *stored_fespace = problem->fespaces.Get("HCurl");

  EXPECT_EQ(stored_fespace->GetVSize(), stored_gf->ParFESpace()->GetVSize());
  ASSERT_TRUE(stored_fespace->GetVSize() > 0);
}
