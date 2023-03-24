#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {
// Specifies output interfaces of a time-domain EM formulation.
class SteadyFormulation : public mfem::Operator {

public:
  SteadyFormulation(){};

  ~SteadyFormulation(){};

  virtual void Init(mfem::Vector &X){};

  virtual void RegisterVariables(){};

  virtual void RegisterOutputFields(mfem::DataCollection *dc_){};

  virtual void WriteOutputFields(mfem::DataCollection *dc_, int it = 0){};

  virtual void WriteConsoleSummary(double t, int it){};

  virtual void InitializeGLVis(){};

  virtual void DisplayToGLVis(){};

  virtual void NotMult(const mfem::Vector& x, mfem::Vector& y){};

  mfem::Array<int> true_offsets;
};
} // namespace hephaestus
