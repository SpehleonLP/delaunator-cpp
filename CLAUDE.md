# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

delaunator-cpp is a header-only C++ library for fast Delaunay triangulation of 2D points, ported from the JavaScript implementation at https://github.com/mapbox/delaunator. The library is designed for performance and can handle millions of points.

## Build System

The project uses CMake with a Makefile wrapper. Key commands:

### Build Commands
- `make release` - Build in release mode (optimized, -O3)
- `make debug` - Build in debug mode (-O0, DEBUG flag)
- `make clean` - Remove build artifacts from cmake-build/
- `make distclean` - Remove all build artifacts including mason dependencies

Build artifacts are placed in `cmake-build/` by default.

### Testing & Benchmarking
- `make test` - Run unit tests (must build first)
- `make bench` - Run benchmarks (must build first)

The test suite uses Catch2 framework, located in `test/delaunator.test.cpp`.

### Code Formatting
- `make format` - Run clang-format on the codebase
- Format settings in `.clang-format`: C++11 style, 4-space indent, no tabs, no column limit

## Project Structure

### Core Library
- `include/delaunator.hpp` - Main header with class declaration and type definitions
- `include/delunator.cpp` - Implementation file with the triangulation algorithm
- The library is header-only in practice (include both files)

### Precision Configuration
The library supports configurable precision via `DelaunatorPrecision` macro:
- `DelaunatorPrecision=64` - Uses `double` and `size_t` (d_fp, d_size types)
- Other values - Uses `float` and `uint32_t`
- The macro **must** be defined before including the header

### Key Data Structures
- `Delaunator` class takes a flat coordinate vector: `{x0, y0, x1, y1, x2, y2, ...}`
- Outputs:
  - `triangles` - flat array of point indices, 3 per triangle
  - `halfedges` - adjacency information between triangles
  - `hull_prev`, `hull_next`, `hull_tri` - convex hull information

### Examples
Located in `examples/`:
- `basic.cpp` - Simple usage demonstration
- `triangulate_geojson.cpp` - Real-world GeoJSON triangulation example
- `utils.hpp` - Helper utilities

### Dependencies
Managed via Mason package manager (mason.cmake):
- Catch2 (v2.4.0) - Testing framework
- RapidJSON (v1.1.0) - JSON parsing for examples
- Google Benchmark (v1.2.0) - Performance benchmarking

## Development Notes

### Compiler Requirements
- C++14 standard required
- Strict warning flags enabled by default (see CMakeLists.txt line 24)
- `-Werror` enabled by default (disable with `WERROR=false make release`)

### Algorithm Details
The implementation uses:
- Spatial hashing for point lookup optimization
- Iterative legalization (non-recursive, using edge stack)
- Kahan-Babuska summation for reduced floating-point error
- Orient and circumradius predicates for geometric tests

### Benchmark Options
CMake options for large-scale benchmarking:
- `-DBENCHMARK_BIG_O=ON` - Calculate algorithmic complexity
- `-DBENCHMARK_100M=ON` - Test with 100M points
- `-DBENCHMARK_10M=ON` - Test with 10M points
