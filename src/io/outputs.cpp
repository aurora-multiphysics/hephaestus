#include "outputs.hpp"

namespace hephaestus {

Outputs::Outputs() {}

Outputs::Outputs(hephaestus::GridFunctions &gridfunctions)
    : _gridfunctions(&gridfunctions), _cycle(0), _use_glvis(false),
      _my_comm(MPI_COMM_WORLD) {
  MPI_Comm_size(_my_comm, &_n_ranks);
  MPI_Comm_rank(_my_comm, &_my_rank);
}
} // namespace hephaestus
