#include "variables.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

TEST(VariablesTest, CheckSetup) {
  mfem::Mesh mesh(
      (std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);
  mfem::ParMesh pmesh = mfem::ParMesh(MPI_COMM_WORLD, mesh);

  hephaestus::Variables variables;
  hephaestus::InputParameters hcurlvarparams;

  hcurlvarparams.SetParam("VariableName", std::string("vector_potential"));
  hcurlvarparams.SetParam("FESpaceName", std::string("HCurl"));
  hcurlvarparams.SetParam("FESpaceType", std::string("ND"));
  hcurlvarparams.SetParam("order", 2);
  hcurlvarparams.SetParam("components", 3);

  variables.AddVariable(hcurlvarparams);
  variables.Init(pmesh);
  mfem::ParGridFunction *stored_gf = variables.gfs.Get("vector_potential");
  mfem::ParFiniteElementSpace *stored_fespace = variables.fespaces.Get("HCurl");

  EXPECT_EQ(stored_fespace->GetVSize(), stored_gf->ParFESpace()->GetVSize());
  ASSERT_TRUE(stored_fespace->GetVSize() > 0);
}
