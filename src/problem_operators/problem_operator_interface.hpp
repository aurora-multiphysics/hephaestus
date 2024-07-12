#pragma once
#include "problem_builder_base.hpp"

namespace hephaestus
{
/// Interface inherited by ProblemOperator and TimeDomainProblemOperator. Removes duplicated code in both classes.
class ProblemOperatorInterface
{
public:
  virtual ~ProblemOperatorInterface() = default;

  virtual void Init() = 0;

  virtual void Update() = 0;

protected:
  virtual void SetTrialVariableNames() = 0;
  virtual void SetTrialVariables() = 0;
  virtual void UpdateOffsets() = 0;

  virtual void UpdateBlockVector(mfem::BlockVector & X) = 0;

  virtual void UpdateOffsetsWithSize(size_t soln_vector_size) = 0;

  /// Returns a reference to the operator's width.
  virtual int & Width() = 0;

  /// Returns a reference to the operator's height.
  virtual int & Height() = 0;
};
}