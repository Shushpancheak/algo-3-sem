#include "ConvertZFuncToPrefixFunc.hpp"
#include "gtest/gtest.h"


namespace {

class ConverterZToPrefixTest : public ::testing::Test {
 protected:
};

TEST_F(ConverterZToPrefixTest, z_to_prefix_1) {
  // abacabad
  std::vector<size_t> z_func = {0, 0, 1, 0, 3, 0, 1, 0};
  std::vector<size_t> p_func = {0, 0, 1, 0, 1, 2, 3, 0};
  ASSERT_EQ(ConvertZFuncToPrefixFunc(z_func), p_func);
  // abbabbdca
  z_func = {0, 0, 0, 3, 0, 0, 0, 0, 1};
  p_func = {0, 0, 0, 1, 2, 3, 0, 0, 1};
  ASSERT_EQ(ConvertZFuncToPrefixFunc(z_func), p_func);
  // abacaadada
  z_func = {0, 0, 1, 0, 1, 1, 0, 1, 0, 1};
  p_func = {0, 0, 1, 0, 1, 1, 0, 1, 0, 1};
  ASSERT_EQ(ConvertZFuncToPrefixFunc(z_func), p_func);
  // abbabbacabbad
  z_func = {0, 0, 0, 4, 0, 0, 1, 0, 4, 0, 0, 1, 0};
  p_func = {0, 0, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0};
  ASSERT_EQ(ConvertZFuncToPrefixFunc(z_func), p_func);
}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}