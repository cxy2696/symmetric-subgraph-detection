// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui

#include <iostream>

#include "test_util.hpp"

void runLdltTests();
void runDetectorTests();
void runGraphTests();

int main() {
    runLdltTests();
    runDetectorTests();
    runGraphTests();

    std::cout << testutil::checks << " checks, " << testutil::failures << " failure(s)\n";
    return testutil::failures == 0 ? 0 : 1;
}
