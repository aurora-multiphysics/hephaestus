#include "hephaestus_steady.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestLEForm : public testing::Test {
protected:
  
  hephaestus::InputParameters test_params() {

    // hephaestus::Subdomain mat1("mat1", 1);
    // mat1.property_map["lambda"] =
    //     new mfem::ConstantCoefficient(1.0);

    // mat1.property_map["shear_mod"] =
    //     new mfem::ConstantCoefficient(1.0);

    // hephaestus::Subdomain mat2("mat2", 2);
    // mat2.property_map["lambda"] =
    //     new mfem::ConstantCoefficient(50.0);

    // mat2.property_map["shear_mod"] =
    //     new mfem::ConstantCoefficient(50.0);

    hephaestus::DomainProperties domain_properties;
        // std::vector<hephaestus::Subdomain>({mat1, mat2}));

    hephaestus::BCMap bc_map;
 
    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./myBeam_3.g")).c_str(),
        1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("AVFormVisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("AVFormParaView");

    hephaestus::Outputs outputs(data_collections);
    hephaestus::Variables variables;
    hephaestus::Postprocessors postprocessors;
    hephaestus::AuxKernels auxkernels;
    hephaestus::Sources sources;

    hephaestus::InputParameters exec_params;
    exec_params.SetParam("VisualisationSteps", int(1));
    exec_params.SetParam("UseGLVis", false);
    hephaestus::SteadyExecutioner *executioner =
        new hephaestus::SteadyExecutioner(exec_params);

    hephaestus::InputParameters params;
    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("Executioner", executioner);
    params.SetParam("Order", 1);
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("DomainProperties", domain_properties);
    params.SetParam("Variables", variables);
    params.SetParam("AuxKernels", auxkernels);
    params.SetParam("Postprocessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);
    params.SetParam("FormulationName", std::string("LinearElasticForm"));

    return params;
  }
};

TEST_F(TestLEForm, CheckRun) {
  hephaestus::InputParameters params(test_params());
  std::cout << "Test begin" << std::endl;
  hephaestus::SteadyExecutioner *executioner(
      params.GetParam<hephaestus::SteadyExecutioner *>("Executioner"));
  executioner->Init(params);
  executioner->Solve();
}
