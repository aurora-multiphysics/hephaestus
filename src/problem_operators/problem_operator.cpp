#include "problem_operator.hpp"

namespace hephaestus
{

void
ProblemOperator::UpdateOperatorWidthAndHeight()
{
  width = height = _true_offsets.Last();
}

} // namespace hephaestus