#include "logging.hpp"

namespace spdlog::sinks
{

template <typename Mutex>
void
MpiSink<Mutex>::sink_it_(const spdlog::details::log_msg & msg)
{
  if (mfem::Mpi::WorldRank() == 0)
  {
    spdlog::memory_buf_t formatted;
    spdlog::sinks::base_sink<Mutex>::formatter_->format(msg, formatted);
    std::cout << fmt::to_string(formatted);
  }
}

template <typename Mutex>
void
MpiSink<Mutex>::flush_()
{
  std::cout << std::flush;
}

}

// Global logger
spdlog::logger logger("Hephaestus Logger", std::make_shared<spdlog::sinks::MpiSink_st>());
