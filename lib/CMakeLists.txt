cmake_minimum_required(VERSION 3.12)
project(mathlib)

set(SOURCE testUtils.h testUtils.cpp
           integration.h integration.cpp)
add_library(${PROJECT_NAME} ${SOURCE})

find_package(fmt CONFIG REQUIRED)
find_package(doctest CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt fmt::fmt-header-only doctest::doctest)

