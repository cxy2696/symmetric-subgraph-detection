cmake_minimum_required(VERSION 3.16)
project(sna VERSION 0.2.0 LANGUAGES CXX)

# Header-only library.
add_library(sna INTERFACE)
add_library(sna::sna ALIAS sna)
target_include_directories(sna INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(sna INTERFACE cxx_std_17)

include(GNUInstallDirs)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

option(SNA_BUILD_EXAMPLES "Build the examples" ON)
option(SNA_BUILD_TESTS "Build the tests" ON)
option(SNA_BUILD_BENCH "Build the benchmark" ON)

if(NOT MSVC)
    add_compile_options(-Wall -Wextra -Wpedantic)
else()
    add_compile_options(/W4)
endif()

if(SNA_BUILD_EXAMPLES)
    add_executable(example_cliques examples/cliques.cpp)
    target_link_libraries(example_cliques PRIVATE sna::sna)

    add_executable(example_extensions examples/extensions.cpp)
    target_link_libraries(example_extensions PRIVATE sna::sna)

    add_executable(example_edge_list examples/edge_list.cpp)
    target_link_libraries(example_edge_list PRIVATE sna::sna)
endif()

if(SNA_BUILD_TESTS)
    enable_testing()
    add_executable(sna_tests
        tests/test_main.cpp
        tests/test_ldlt.cpp
        tests/test_detector.cpp
        tests/test_graph.cpp
    )
    target_link_libraries(sna_tests PRIVATE sna::sna)
    add_test(NAME sna_tests COMMAND sna_tests)
endif()

if(SNA_BUILD_BENCH)
    add_executable(bench_scale bench/bench.cpp)
    target_link_libraries(bench_scale PRIVATE sna::sna)
endif()
