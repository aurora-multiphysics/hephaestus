#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.hpp"

namespace hephaestus {
// Specifies output interfaces of a time-domain EM formulation.
class TransientFormulation : public mfem::TimeDependentOperator {

public:
  TransientFormulation(){};

  ~TransientFormulation(){};

  virtual void Init(mfem::Vector &X){};

  virtual void RegisterOutputFields(mfem::DataCollection *dc_){};

  virtual void WriteOutputFields(mfem::DataCollection *dc_, int it = 0){};

  virtual void WriteConsoleSummary(double t, int it){};

  virtual void InitializeGLVis(){};

  virtual void DisplayToGLVis(){};

  mfem::Array<int> true_offsets;
};
} // namespace hephaestus
