#include "inputs.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

TEST(InputParametersTest, CheckData) {
  hephaestus::InputParameters params;
  int example_int = 5;
  params.SetParam("IntegerParam", example_int);

  std::string example_string("ExampleString");
  params.SetParam("StringParam", example_string);

  mfem::Array<int> example_array({1, 2, 3});
  params.SetParam("ArrayParam", example_array);

  EXPECT_EQ(params.GetParam<int>("IntegerParam"), example_int);

  EXPECT_EQ(params.GetParam<std::string>("StringParam"), example_string);

  mfem::Array<int> stored_array =
      params.GetParam<mfem::Array<int>>("ArrayParam");
  for (int i = 0; i < example_array.Size(); ++i) {
    EXPECT_EQ(example_array[i], stored_array[i])
        << "Input and stored arrays differ at index " << i;
  }
}
