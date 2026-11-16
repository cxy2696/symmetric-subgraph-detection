// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// Minimal assertion helpers. The project intentionally has no test framework
// dependency: the test binary is enough for `ctest`.

#ifndef SNA_TEST_UTIL_HPP
#define SNA_TEST_UTIL_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace testutil {

inline int checks = 0;
inline int failures = 0;

inline void record(bool ok, const char* expression, const char* file, int line) {
    ++checks;
    if (ok) return;
    ++failures;
    std::cerr << file << ":" << line << ": check failed: " << expression << "\n";
}

template <class T>
bool nearlyEqual(T a, T b, T tolerance) {
    const T scale = std::max(T(1), std::max(std::abs(a), std::abs(b)));
    return std::abs(a - b) <= tolerance * scale;
}

inline std::vector<std::vector<int>> canonical(std::vector<std::vector<int>> groups) {
    for (std::vector<int>& group : groups) std::sort(group.begin(), group.end());
    std::sort(groups.begin(), groups.end());
    return groups;
}

inline std::string format(const std::vector<std::vector<int>>& groups) {
    std::string out = "[";
    for (std::size_t i = 0; i < groups.size(); ++i) {
        if (i != 0) out += ", ";
        out += "{";
        for (std::size_t j = 0; j < groups[i].size(); ++j) {
            if (j != 0) out += ",";
            out += std::to_string(groups[i][j]);
        }
        out += "}";
    }
    return out + "]";
}

inline void checkGroups(const std::vector<std::vector<int>>& got,
                        const std::vector<std::vector<int>>& expected,
                        const char* file, int line) {
    ++checks;
    const std::vector<std::vector<int>> a = canonical(got);
    const std::vector<std::vector<int>> b = canonical(expected);
    if (a == b) return;
    ++failures;
    std::cerr << file << ":" << line << ": group lists differ\n"
              << "  got      " << format(a) << "\n"
              << "  expected " << format(b) << "\n";
}

inline void checkThrows(bool threw, const char* expression, const char* file, int line) {
    ++checks;
    if (threw) return;
    ++failures;
    std::cerr << file << ":" << line << ": expected exception: " << expression << "\n";
}

}  // namespace testutil

#define CHECK(condition) ::testutil::record((condition), #condition, __FILE__, __LINE__)
#define CHECK_NEAR(a, b, tolerance) \
    ::testutil::record(::testutil::nearlyEqual((a), (b), (tolerance)), #a " ~= " #b, __FILE__, __LINE__)
#define CHECK_GROUPS(got, expected) ::testutil::checkGroups((got), (expected), __FILE__, __LINE__)
#define CHECK_THROWS(statement)                                                    \
    do {                                                                           \
        bool threw_ = false;                                                       \
        try {                                                                      \
            statement;                                                             \
        } catch (...) {                                                            \
            threw_ = true;                                                         \
        }                                                                          \
        ::testutil::checkThrows(threw_, #statement, __FILE__, __LINE__);            \
    } while (false)

#endif  // SNA_TEST_UTIL_HPP
