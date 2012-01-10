#include <gtest/gtest.h>

#include "unittests_common.hh"
#include "unittests_iterators.hh"

int main(int _argc, char** _argv) {

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();
}
