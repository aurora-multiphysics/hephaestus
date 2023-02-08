#include "variables.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

TEST(VariablesTest, CheckSetup) {
  mfem::Mesh mesh(
      (std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);
  mfem::ParMesh pmesh = mfem::ParMesh(MPI_COMM_WORLD, mesh);

  hephaestus::FESpaces fespaces;
  hephaestus::GridFunctions gridfunctions;

  hephaestus::InputParameters hcurlfespaceparams;
  hcurlfespaceparams.SetParam("FESpaceName", std::string("HCurl"));
  hcurlfespaceparams.SetParam("FESpaceType", std::string("ND"));
  hcurlfespaceparams.SetParam("order", 2);
  hcurlfespaceparams.SetParam("components", 3);
  fespaces.StoreInput(hcurlfespaceparams);

  hephaestus::InputParameters vectorpotentialparams;
  vectorpotentialparams.SetParam("VariableName",
                                 std::string("vector_potential"));
  vectorpotentialparams.SetParam("FESpaceName", std::string("HCurl"));
  gridfunctions.StoreInput(vectorpotentialparams);

  fespaces.Init(pmesh);
  gridfunctions.Init(pmesh, fespaces);

  mfem::ParGridFunction *stored_gf = gridfunctions.Get("vector_potential");
  mfem::ParFiniteElementSpace *stored_fespace = fespaces.Get("HCurl");

  EXPECT_EQ(stored_fespace->GetVSize(), stored_gf->ParFESpace()->GetVSize());
  ASSERT_TRUE(stored_fespace->GetVSize() > 0);
}
