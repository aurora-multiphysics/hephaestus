#include "steady_formulation.hpp"

namespace hephaestus {

void FrequencyDomainOperator::SetVariables() {
  // state_var_names.push_back("electric_field_real");
  // state_var_names.push_back("electric_field_imag");

  local_test_vars = populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
      _variables, state_var_names);
  e_ = new mfem::ParComplexGridFunction(local_test_vars.at(0)->ParFESpace());
  *e_ = std::complex(0.0, 0.0);

  // Set operator size and block structure
  block_trueOffsets.SetSize(local_test_vars.size() + 1);
  block_trueOffsets[0] = 0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    block_trueOffsets[ind + 1] =
        local_test_vars.at(ind)->ParFESpace()->TrueVSize();
  }
  block_trueOffsets.PartialSum();

  true_offsets.SetSize(local_test_vars.size() + 1);
  true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    true_offsets[ind + 1] = local_test_vars.at(ind)->ParFESpace()->GetVSize();
  }
  true_offsets.PartialSum();

  this->height = true_offsets[local_test_vars.size()];
  this->width = true_offsets[local_test_vars.size()];
  trueX.Update(block_trueOffsets);
  trueRhs.Update(block_trueOffsets);

  // Populate vector of active auxiliary variables
  active_aux_var_names.resize(0);
  for (auto &aux_var_name : aux_var_names) {
    if (_variables.Has(aux_var_name)) {
      active_aux_var_names.push_back(aux_var_name);
    }
  }
};

void FrequencyDomainOperator::Init(mfem::Vector &X) {
  // Define material property coefficients
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    local_test_vars.at(ind)->MakeRef(local_test_vars.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
  }

  // muInvCoef_ =
  // _domain_properties.scalar_property_map["magnetic_reluctivity"]; massCoef_ =
  // _domain_properties.scalar_property_map["hertz_mass"]; lossCoef_ =
  // _domain_properties.scalar_property_map["hertz_loss"];
}

void FrequencyDomainOperator::Solve(mfem::Vector &X) {
  // mfem::OperatorHandle A1;
  // mfem::Vector E, RHS;
  // mfem::OperatorHandle PCOp;

  // mfem::Vector zeroVec(3);
  // zeroVec = 0.0;
  // mfem::VectorConstantCoefficient zeroCoef(zeroVec);

  // a1_ = new mfem::ParSesquilinearForm(e_->ParFESpace(), conv_);
  // a1_->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*muInvCoef_), NULL);
  // a1_->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*massCoef_),
  // NULL); if (lossCoef_) {
  //   a1_->AddDomainIntegrator(NULL,
  //                            new mfem::VectorFEMassIntegrator(*lossCoef_));
  // }

  // // Volume Current Density
  // mfem::Vector j(3);
  // j = 0.0;
  // jrCoef_ = new mfem::VectorConstantCoefficient(j);
  // jiCoef_ = new mfem::VectorConstantCoefficient(j);

  // jd_ = new mfem::ParComplexLinearForm(e_->ParFESpace(), conv_);
  // jd_->AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*jrCoef_),
  //                          new mfem::VectorFEDomainLFIntegrator(*jiCoef_));

  // _bc_map.applyEssentialBCs(std::string("electric_field"), ess_bdr_tdofs_,
  // *e_,
  //                           pmesh_);
  // _bc_map.applyIntegratedBCs(std::string("electric_field"), *jd_, pmesh_);
  // _bc_map.applyIntegratedBCs(std::string("electric_field"), *a1_, pmesh_);

  // a1_->Assemble();
  // a1_->Finalize();

  // jd_->Assemble();

  // a1_->FormLinearSystem(ess_bdr_tdofs_, *e_, *jd_, A1, E, RHS);

  // mfem::ComplexHypreParMatrix *A1Z = A1.As<mfem::ComplexHypreParMatrix>();
  // mfem::HypreParMatrix *A1C = A1Z->GetSystemMatrix();
  // mfem::SuperLURowLocMatrix A_SuperLU(*A1C);
  // mfem::SuperLUSolver solver(MPI_COMM_WORLD);
  // solver.SetOperator(A_SuperLU);
  // solver.Mult(RHS, E);
  // delete A1C;

  // e_->Distribute(E);

  // *_variables.Get(state_var_names.at(0)) = e_->real();
  // *_variables.Get(state_var_names.at(1)) = e_->imag();
}

SteadyFormulation::SteadyFormulation() : fd_operator(NULL), oneCoef(1.0){};

} // namespace hephaestus
