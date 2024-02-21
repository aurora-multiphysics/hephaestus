#include "mfem.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/base_sink.h"

namespace spdlog::sinks
{

template <typename Mutex>
class MpiSink : public spdlog::sinks::base_sink<Mutex>
{

protected:
  void sink_it_(const spdlog::details::log_msg & msg) override
  {

    // log_msg is a struct containing the log entry info like level, timestamp, thread id etc.
    // msg.raw contains pre formatted log

    // If needed (very likely but not mandatory), the sink formats the message before sending it to
    // its final destination:
    if (mfem::Mpi::WorldRank() == 0)
    {

      spdlog::memory_buf_t formatted;
      spdlog::sinks::base_sink<Mutex>::formatter_->format(msg, formatted);
      std::cout << fmt::to_string(formatted);
    }
  }

  void flush_() override { std::cout << std::flush; }
};

#include "spdlog/details/null_mutex.h"
#include <mutex>
using MpiSink_mt = MpiSink<std::mutex>;
using MpiSink_st = MpiSink<spdlog::details::null_mutex>;

} // namespace spdlog::sinks