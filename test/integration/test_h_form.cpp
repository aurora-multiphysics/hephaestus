#include "hephaestus_h_form.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

TEST(TestHForm, CheckRun) {
  std::vector<char *> argv;
  h_form_solve(0, argv.data());
}
