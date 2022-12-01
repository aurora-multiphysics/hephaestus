// // Based on an H form MMS test provided by Joseph Dean

// #include "auxkernels.hpp"
// #include "executioner.hpp"

// #include "hephaestus_transient.hpp"
// #include "postprocessors.hpp"

// #include <gtest/gtest.h>

// extern const char *DATA_DIR;

// class TestSolenoid : public testing::Test {
// protected:
//   static double estimate_convergence_rate(HYPRE_BigInt n_i, HYPRE_BigInt
//   n_imo,
//                                           double error_i, double error_imo,
//                                           int dim) {
//     return std::log(error_i / error_imo) /
//            std::log(std::pow(n_imo / static_cast<double>(n_i), 1.0 / dim));
//   }

//   static double potential_ground(const mfem::Vector &x, double t) {
//     return -x(0);
//   }

//   static void adot_bc(const mfem::Vector &x, double t, mfem::Vector &H) {
//     H(0) = 1 + sin(x(1) * M_PI) * sin(x(2) * M_PI);
//     H(1) = 0;
//     H(2) = 0;
//   }

//   static void A_exact_expr(const mfem::Vector &x, double t,
//                            mfem::Vector &A_exact) {
//     A_exact(0) = (1 + sin(x(1) * M_PI) * sin(x(2) * M_PI)) * t;
//     A_exact(1) = 0;
//     A_exact(2) = 0;
//   }
//   static double mu_expr(const mfem::Vector &x) {
//     double variation_scale = 0.0;
//     double mu =
//         1.0 / (1.0 + variation_scale * cos(M_PI * x(0)) * cos(M_PI * x(1)));
//     return mu;
//   }
//   static void source_current(const mfem::Vector &x, double t, mfem::Vector
//   &J) {
//     J(0) = -(x(1) / sqrt(x(0) * x(0) + x(1) * x(1))) * t;
//     J(1) = (x(0) / sqrt(x(0) * x(0) + x(1) * x(1))) * t;
//     J(2) = 0.0;
//   }

//   hephaestus::InputParameters test_params() {
//     hephaestus::Subdomain coil("coil", 1);
//     coil.property_map["electrical_conductivity"] =
//         new mfem::ConstantCoefficient(1.0);
//     hephaestus::Subdomain ring("ring", 3);
//     ring.property_map["electrical_conductivity"] =
//         new mfem::ConstantCoefficient(2.273e7);
//     hephaestus::Subdomain air("air", 2);
//     air.property_map["electrical_conductivity"] =
//         new mfem::ConstantCoefficient(1.0);

//     hephaestus::DomainProperties domain_properties(
//         std::vector<hephaestus::Subdomain>({coil, air, ring}));

//     domain_properties.scalar_property_map["magnetic_permeability"] =
//         new mfem::ConstantCoefficient(M_PI * 4.0e-7);

//     hephaestus::BCMap bc_map;
//     // mfem::VectorFunctionCoefficient *adotVecCoef =
//     //     new mfem::VectorFunctionCoefficient(3, adot_bc);
//     // bc_map["tangential_dAdt"] = new hephaestus::VectorFunctionDirichletBC(
//     //     std::string("magnetic_vector_potential"), mfem::Array<int>({1, 2,
//     //     3}), adotVecCoef);
//     // domain_properties.vector_property_map["surface_tangential_dAdt"] =
//     //     adotVecCoef;
//     // domain_properties.scalar_property_map["electrical_conductivity"] =
//     //     new mfem::ConstantCoefficient(1.0);

//     // mfem::Array<int> ground_terminal(1);
//     // ground_terminal[0] = 1;
//     // mfem::FunctionCoefficient *ground_coeff =
//     //     new mfem::FunctionCoefficient(potential_ground);
//     // bc_map["ground_potential"] = new hephaestus::FunctionDirichletBC(
//     //     std::string("electric_potential"), mfem::Array<int>({1, 2, 3}),
//     //     ground_coeff);
//     // domain_properties.scalar_property_map["ground_potential"] =
//     ground_coeff;

//     mfem::VectorFunctionCoefficient *JSrcCoef =
//         new mfem::VectorFunctionCoefficient(3, source_current);

//     mfem::Array<int> source_coil(1);
//     source_coil[0] = 1;
//     mfem::Array<mfem::VectorCoefficient *> sources(1);
//     sources[0] = JSrcCoef;

//     mfem::PWVectorCoefficient *JRestricted =
//         new mfem::PWVectorCoefficient(3, source_coil, sources);

//     // mfem::VectorRestrictedCoefficient *JRestricted =
//     //     new mfem::VectorRestrictedCoefficient(*JSrcCoef, source_coil);
//     domain_properties.vector_property_map["source"] = JRestricted;

//     mfem::VectorFunctionCoefficient *A_exact =
//         new mfem::VectorFunctionCoefficient(3, A_exact_expr);

//     domain_properties.vector_property_map["a_exact_coeff"] = A_exact;

//     mfem::Mesh mesh(
//         (std::string(DATA_DIR) + std::string("./solenoid.g")).c_str(), 1, 1);

//     std::map<std::string, mfem::DataCollection *> data_collections;
//     data_collections["VisItDataCollection"] =
//         new mfem::VisItDataCollection("AVFormVisIt");
//     data_collections["ParaViewDataCollection"] =
//         new mfem::ParaViewDataCollection("AVFormParaView");

//     hephaestus::Outputs outputs(data_collections);

//     hephaestus::InputParameters hcurlvarparams;
//     hcurlvarparams.SetParam("VariableName",
//                             std::string("analytic_vector_potential"));
//     hcurlvarparams.SetParam("FESpaceName", std::string("HCurl"));
//     hcurlvarparams.SetParam("FESpaceType", std::string("ND"));
//     hcurlvarparams.SetParam("order", 4);
//     hcurlvarparams.SetParam("components", 3);
//     hephaestus::Variables variables;
//     variables.AddVariable(hcurlvarparams);

//     hephaestus::InputParameters l2errpostprocparams;
//     l2errpostprocparams.SetParam("VariableName",
//                                  std::string("magnetic_vector_potential"));
//     l2errpostprocparams.SetParam("VectorCoefficientName",
//                                  std::string("a_exact_coeff"));
//     hephaestus::Postprocessors postprocessors;
//     postprocessors.Register(
//         "L2ErrorPostprocessor",
//         new hephaestus::L2ErrorVectorPostprocessor(l2errpostprocparams),
//         true);

//     hephaestus::InputParameters vectorcoeffauxparams;
//     vectorcoeffauxparams.SetParam("VariableName",
//                                   std::string("analytic_vector_potential"));
//     vectorcoeffauxparams.SetParam("VectorCoefficientName",
//                                   std::string("a_exact_coeff"));

//     hephaestus::AuxKernels auxkernels;
//     auxkernels.Register(
//         "VectorCoefficientAuxKernel",
//         new hephaestus::VectorCoefficientAuxKernel(vectorcoeffauxparams),
//         true);

//     hephaestus::InputParameters exec_params;
//     exec_params.SetParam("TimeStep", float(0.05));
//     exec_params.SetParam("StartTime", float(0.00));
//     exec_params.SetParam("EndTime", float(0.05));
//     hephaestus::TransientExecutioner *executioner =
//         new hephaestus::TransientExecutioner(exec_params);

//     hephaestus::InputParameters params;
//     params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
//     params.SetParam("Executioner", executioner);
//     params.SetParam("Order", 2);
//     params.SetParam("BoundaryConditions", bc_map);
//     params.SetParam("DomainProperties", domain_properties);
//     params.SetParam("Variables", variables);
//     params.SetParam("AuxKernels", auxkernels);
//     params.SetParam("Postprocessors", postprocessors);
//     params.SetParam("Outputs", outputs);
//     params.SetParam("FormulationName", std::string("AVForm"));

//     return params;
//   }
// };

// TEST_F(TestSolenoid, CheckRun) {
//   hephaestus::InputParameters params(test_params());

//   hephaestus::TransientExecutioner *executioner(
//       params.GetParam<hephaestus::TransientExecutioner *>("Executioner"));
//   executioner->Init(params);
//   executioner->Solve();

//   hephaestus::L2ErrorVectorPostprocessor l2errpostprocessor =
//       *(dynamic_cast<hephaestus::L2ErrorVectorPostprocessor *>(
//           params.GetParam<hephaestus::Postprocessors>("Postprocessors")
//               .Get("L2ErrorPostprocessor")));

//   double r;
//   for (std::size_t i = 1; i < l2errpostprocessor.ndofs.Size(); ++i) {
//     r = estimate_convergence_rate(
//         l2errpostprocessor.ndofs[i], l2errpostprocessor.ndofs[i - 1],
//         l2errpostprocessor.l2_errs[i], l2errpostprocessor.l2_errs[i - 1], 3);
//     std::cout << r << std::endl;
//     ASSERT_TRUE(r > params.GetParam<int>("Order"));
//     ASSERT_TRUE(r < params.GetParam<int>("Order") + 1.0);
//   }
// }
