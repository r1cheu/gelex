#!/usr/bin/env bash
set -euo pipefail

INPUT=$(cat)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty')

if [[ -z "$FILE_PATH" ]]; then
    exit 0
fi

case "$FILE_PATH" in
    *.cpp|*.h|*.hpp|*.cc|*.cxx)
        if command -v clang-format &>/dev/null && [[ -f "$FILE_PATH" ]]; then
            clang-format -i "$FILE_PATH"
        fi
        ;;
esac
