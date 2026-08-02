#!/usr/bin/env bash
# Copyright 2026 RuLei Chen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
