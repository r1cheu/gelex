#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <version>"
    echo "Example: $0 0.9.0"
    exit 1
fi

VERSION="$1"

if [[ ! "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "Error: Version must be in format X.Y.Z (e.g., 0.9.0)"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

echo "Updating version to $VERSION..."

sed -i "s/VERSION [0-9]\+\.[0-9]\+\.[0-9]\+/VERSION $VERSION/" CMakeLists.txt
sed -i "s/version: \"[0-9]\+\.[0-9]\+\.[0-9]\+\"/version: \"$VERSION\"/" recipe.yaml

echo "Creating git commit and tag..."

git add CMakeLists.txt recipe.yaml
git commit -m "chore: bump version to $VERSION"
git tag -a "v$VERSION" -m "Release v$VERSION"

echo ""
echo "Done! Version updated to $VERSION"
echo ""
echo "To push the release:"
echo "  git push && git push --tags"
