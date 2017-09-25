#include "gtest/gtest.h"
#include <DeformetricaConfig.h>

std::vector<std::string> args;

int main(int argc, char **argv) {
    args = std::vector<std::string>(argv, argv + argc);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}