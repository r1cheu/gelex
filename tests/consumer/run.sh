#!/usr/bin/env bash
# Copyright 2026 RuLei Chen
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
install_prefix="${project_root}/build/ci-install"
consumer_build_dir="${project_root}/build/ci-consumer"

cmake -E remove_directory "${install_prefix}"
cmake -E remove_directory "${consumer_build_dir}"
cmake --install "${project_root}/build/ci" --prefix "${install_prefix}"
cmake \
    -S "${project_root}/tests/consumer" \
    -B "${consumer_build_dir}" \
    -G Ninja \
    -DCMAKE_PREFIX_PATH="${install_prefix}"
cmake --build "${consumer_build_dir}"
"${consumer_build_dir}/gelex_consumer"
