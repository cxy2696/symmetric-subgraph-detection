cmake_minimum_required(VERSION 3.16)
project(sna VERSION 0.3.0 LANGUAGES CXX)

include(GNUInstallDirs)

# Header-only library.
add_library(sna INTERFACE)
add_library(sna::sna ALIAS sna)
target_include_directories(sna INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(sna INTERFACE cxx_std_17)

include(CMakePackageConfigHelpers)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS sna EXPORT snaTargets)
install(EXPORT snaTargets
        FILE snaTargets.cmake
        NAMESPACE sna::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sna)
configure_package_config_file(
    cmake/snaConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/snaConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sna)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/snaConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion)
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/snaConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/snaConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sna)

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
