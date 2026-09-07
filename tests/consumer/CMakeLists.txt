# cmake-format: off
# Copyright 2026 RuLei Chen
# SPDX-License-Identifier: Apache-2.0
# cmake-format: on

cmake_minimum_required(VERSION 3.20)

project(gelex_consumer LANGUAGES CXX)

find_package(gelex CONFIG REQUIRED)

add_executable(gelex_consumer main.cpp)
target_link_libraries(gelex_consumer PRIVATE gelex::gelex)
