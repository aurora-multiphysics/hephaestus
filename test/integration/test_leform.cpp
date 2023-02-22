#include "hephaestus_steady.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestLEForm : public testing::Test {
protected:
  
  hephaestus::InputParameters test_params() {

    hephaestus::Subdomain mat1("mat1", 1);
    mat1.property_map["lambda"] =
        new mfem::ConstantCoefficient(50.0);

    mat1.property_map["shear_modulus"] =
        new mfem::ConstantCoefficient(50.0);

    hephaestus::Subdomain mat2("mat2", 2);
    mat2.property_map["lambda"] =
        new mfem::ConstantCoefficient(1.0);

    mat2.property_map["shear_modulus"] =
        new mfem::ConstantCoefficient(1.0);

    hephaestus::DomainProperties domain_properties(
        std::vector<hephaestus::Subdomain>({mat1, mat2}));

    
    mfem::VectorArrayCoefficient f(3);
    // Set all components to 0
    for (int i = 0; i < 2; i++)
    {
        f.Set(i, new mfem::ConstantCoefficient(0.0));
    }
    {
        f.Set(2, new mfem::ConstantCoefficient(-1.0e-2));
    }       

    domain_properties.vector_property_map["push_force"] = new mfem::VectorArrayCoefficient(3);
    dynamic_cast<mfem::VectorArrayCoefficient *>(domain_properties.vector_property_map["push_force"])->Set(0, new mfem::ConstantCoefficient(0.0));
    dynamic_cast<mfem::VectorArrayCoefficient *>(domain_properties.vector_property_map["push_force"])->Set(1, new mfem::ConstantCoefficient(0.0));
    dynamic_cast<mfem::VectorArrayCoefficient *>(domain_properties.vector_property_map["push_force"])->Set(2, new mfem::ConstantCoefficient(-1.0e-2));


    hephaestus::BCMap bc_map;

    mfem::Array<int> bdr_attrs(1);
    bdr_attrs[0] = 2;

    // mfem::VectorBoundaryLFIntegrator *integ = new mfem::VectorBoundaryLFIntegrator(f);
    mfem::VectorBoundaryLFIntegrator *integ = new mfem::VectorBoundaryLFIntegrator(*domain_properties.vector_property_map["push_force"]);
    hephaestus::IntegratedBC *bc = new hephaestus::IntegratedBC(std::string("displacement"), bdr_attrs, integ);
    bc_map["linear_elastic_force_bc"] = bc;

    // mfem::Array<int> wgi_in_attr(1);
    // wgi_in_attr[0] = 2;
    // hephaestus::IntegratedBC waveguide_in(std::string("displacement"),
                                            // wgi_in_attr);
    // mfem::VectorFunctionCoefficient UReal(3, f);
    // mfem::VectorFunctionCoefficient UImag(pmesh.SpaceDimension(), RWTE10_imag);
    // waveguide_in.lfi_re = new mfem::VectorBoundaryLFIntegrator(f);

    // bc_map["linear_elastic_force_bc"] = new hephaestus::IntegratedBC(waveguide_in);



 
    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./myBeam_4.g")).c_str(),
        1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    // data_collections["VisItDataCollection"] =
    //     new mfem::VisItDataCollection("LEFormVisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("LEFormParaView");

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
