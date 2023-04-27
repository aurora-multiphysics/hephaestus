#include "auxkernels.hpp"
#include "steady_executioner.hpp"

#include "hephaestus_transient.hpp"
#include "postprocessors.hpp"

#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestComplexTeam7 : public testing::Test {
protected:
  static void source_current(const mfem::Vector &xv, mfem::Vector &J) {
    double x0(194e-3);  // Coil centre x coordinate
    double y0(100e-3);  // Coil centre y coordinate
    double a(50e-3);    // Coil thickness
    double I0(2742);    // Coil current in Ampere-turns
    double S(2.5e-3);   // Coil cross sectional area
    double freq(200.0); // Frequency in Hz

    double x = xv(0);
    double y = xv(1);

    // Current density magnitude
    double Jmag = (I0 / S);

    // Calculate x component of current density unit vector
    if (abs(x - x0) < a) {
      J(0) = -(y - y0) / abs(y - y0);
    } else if (abs(y - y0) < a) {
      J(0) = 0.0;
    } else {
      J(0) = -(y - (y0 + a * ((y - y0) / abs(y - y0)))) /
             hypot(x - (x0 + a * ((x - x0) / abs(x - x0))),
                   y - (y0 + a * ((y - y0) / abs(y - y0))));
    }

    // Calculate y component of current density unit vector
    if (abs(y - y0) < a) {
      J(1) = (x - x0) / abs(x - x0);
    } else if (abs(x - x0) < a) {
      J(1) = 0.0;
    } else {
      J(1) = (x - (x0 + a * ((x - x0) / abs(x - x0)))) /
             hypot(x - (x0 + a * ((x - x0) / abs(x - x0))),
                   y - (y0 + a * ((y - y0) / abs(y - y0))));
    }

    // Calculate z component of current density unit vector
    J(2) = 0.0;

    // Scale by current density magnitude
    J *= Jmag;
  }

  hephaestus::InputParameters test_params() {
    hephaestus::Subdomain air("air", 1);
    air.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain plate("plate", 2);
    plate.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(3.526e7);
    hephaestus::Subdomain coil1("coil1", 3);
    coil1.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain coil2("coil2", 4);
    coil2.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain coil3("coil3", 5);
    coil3.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain coil4("coil4", 6);
    coil4.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);

    hephaestus::DomainProperties domain_properties(
        std::vector<hephaestus::Subdomain>(
            {air, plate, coil1, coil2, coil3, coil4}));

    domain_properties.scalar_property_map["frequency"] =
        new mfem::ConstantCoefficient(200.0);
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::ConstantCoefficient(M_PI * 4.0e-7);
    hephaestus::BCMap bc_map;

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./team7_small.g")).c_str(), 1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("Team7VisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("Team7ParaView");
    hephaestus::Outputs outputs(data_collections);

    hephaestus::InputParameters hcurlfespaceparams;
    hcurlfespaceparams.SetParam("FESpaceName", std::string("HCurl"));
    hcurlfespaceparams.SetParam("FESpaceType", std::string("ND"));
    hcurlfespaceparams.SetParam("order", 1);
    hcurlfespaceparams.SetParam("components", 3);
    hephaestus::InputParameters h1fespaceparams;
    h1fespaceparams.SetParam("FESpaceName", std::string("H1"));
    h1fespaceparams.SetParam("FESpaceType", std::string("H1"));
    h1fespaceparams.SetParam("order", 1);
    h1fespaceparams.SetParam("components", 3);
    hephaestus::FESpaces fespaces;
    fespaces.StoreInput(hcurlfespaceparams);
    fespaces.StoreInput(h1fespaceparams);

    hephaestus::InputParameters afieldparams_re;
    afieldparams_re.SetParam("VariableName",
                             std::string("magnetic_vector_potential_real"));
    afieldparams_re.SetParam("FESpaceName", std::string("HCurl"));
    hephaestus::InputParameters afieldparams_im;
    afieldparams_im.SetParam("VariableName",
                             std::string("magnetic_vector_potential_imag"));
    afieldparams_im.SetParam("FESpaceName", std::string("HCurl"));
    hephaestus::InputParameters bfieldparams_re;
    bfieldparams_re.SetParam("VariableName",
                             std::string("magnetic_flux_density_real"));
    bfieldparams_re.SetParam("FESpaceName", std::string("HCurl"));
    hephaestus::InputParameters bfieldparams_im;
    bfieldparams_im.SetParam("VariableName",
                             std::string("magnetic_flux_density_imag"));
    bfieldparams_im.SetParam("FESpaceName", std::string("HCurl"));
    hephaestus::GridFunctions gridfunctions;
    gridfunctions.StoreInput(afieldparams_re);
    gridfunctions.StoreInput(afieldparams_im);
    gridfunctions.StoreInput(bfieldparams_re);
    gridfunctions.StoreInput(bfieldparams_im);

    hephaestus::Postprocessors postprocessors;
    hephaestus::AuxKernels auxkernels;

    hephaestus::Sources sources;
    mfem::VectorFunctionCoefficient *JSrcCoef =
        new mfem::VectorFunctionCoefficient(3, source_current);
    mfem::Array<mfem::VectorCoefficient *> sourcecoefs(4);
    sourcecoefs[0] = JSrcCoef;
    sourcecoefs[1] = JSrcCoef;
    sourcecoefs[2] = JSrcCoef;
    sourcecoefs[3] = JSrcCoef;
    mfem::Array<int> coilsegments(4);
    coilsegments[0] = 3;
    coilsegments[1] = 4;
    coilsegments[2] = 5;
    coilsegments[3] = 6;
    mfem::PWVectorCoefficient *JSrcRestricted =
        new mfem::PWVectorCoefficient(3, coilsegments, sourcecoefs);
    domain_properties.vector_property_map["source"] = JSrcRestricted;

    hephaestus::InputParameters div_free_source_params;
    div_free_source_params.SetParam("SourceName", std::string("source"));
    div_free_source_params.SetParam("HCurlFESpaceName", std::string("HCurl"));
    div_free_source_params.SetParam("H1FESpaceName", std::string("H1"));
    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-12));
    current_solver_options.SetParam("MaxIter", (unsigned int)200);
    current_solver_options.SetParam("PrintLevel", 1);
    div_free_source_params.SetParam("SolverOptions", current_solver_options);
    sources.Register(
        "source",
        new hephaestus::DivFreeVolumetricSource(div_free_source_params), true);

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);
    solver_options.SetParam("PrintLevel", 0);

    hephaestus::SteadyFormulation *formulation =
        new hephaestus::ComplexAFormulation();

    hephaestus::InputParameters params;
    params.SetParam("UseGLVis", true);

    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("Order", 2);
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("DomainProperties", domain_properties);
    params.SetParam("FESpaces", fespaces);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("AuxKernels", auxkernels);
    params.SetParam("Postprocessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);
    params.SetParam("Formulation", formulation);
    params.SetParam("SolverOptions", solver_options);
    std::cout << "Created params ";
    return params;
  }
};

TEST_F(TestComplexTeam7, CheckRun) {
  hephaestus::InputParameters params(test_params());

  hephaestus::SteadyExecutioner *executioner =
      new hephaestus::SteadyExecutioner(params);
  std::cout << "Created exec ";
  executioner->Init();
  executioner->Solve();
}
