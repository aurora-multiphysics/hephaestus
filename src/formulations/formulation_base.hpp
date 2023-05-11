#pragma once
#include "../common/pfem_extras.hpp"
#include "auxkernels.hpp"
#include "equation_system.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

// Specifies output of an EM formulation.
class Formulation {

public:
  mfem::ConstantCoefficient oneCoef{1.0};
  Formulation(){};
  virtual hephaestus::EquationSystem *CreateEquationSystem() const = 0;
  virtual mfem::Operator *
  CreateOperator(mfem::ParMesh &pmesh,
                 mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
                 mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                 hephaestus::BCMap &bc_map,
                 hephaestus::DomainProperties &domain_properties,
                 hephaestus::Sources &sources,
                 hephaestus::InputParameters &solver_options) const = 0;

  virtual void RegisterMissingVariables(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables){};

  virtual void
  RegisterAuxKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxKernels &auxkernels){};

  virtual void
  RegisterCoefficients(hephaestus::DomainProperties &domain_properties){};
};
} // namespace hephaestus
