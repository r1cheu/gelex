#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


INCLUDE_RE = re.compile(r"^\s*#\s*include\s*[<\"]([^\">]+)[\">]")
REL_SRC_INCLUDE_RE = re.compile(r"^(?:\.\./)+src/|^src/")


@dataclass
class Violation:
    rule_id: str
    file: str
    line: int
    include: str
    message: str


def as_posix(path: Path) -> str:
    return path.as_posix()


def rel_to_root(path: Path, root: Path) -> str:
    return as_posix(path.relative_to(root))


def collect_files(root: Path) -> list[Path]:
    scan_roots = [
        root / "include",
        root / "src",
        root / "apps",
        root / "tests",
        root / "benchmark",
    ]
    suffixes = {".h", ".hpp", ".hh", ".hxx", ".c", ".cc", ".cpp", ".cxx"}

    files: list[Path] = []
    for scan_root in scan_roots:
        if not scan_root.exists():
            continue
        for path in scan_root.rglob("*"):
            if path.is_file() and path.suffix.lower() in suffixes:
                files.append(path)
    return files


def iter_includes(path: Path) -> list[tuple[int, str]]:
    entries: list[tuple[int, str]] = []
    text = path.read_text(encoding="utf-8", errors="ignore")
    for line_no, line in enumerate(text.splitlines(), start=1):
        match = INCLUDE_RE.match(line)
        if match:
            entries.append((line_no, match.group(1).strip()))
    return entries


def starts_with(rel_file: str, prefix: str) -> bool:
    return rel_file == prefix or rel_file.startswith(prefix + "/")


def is_rel_src_include(include_target: str) -> bool:
    return bool(REL_SRC_INCLUDE_RE.match(include_target))


def check_boundaries(root: Path) -> list[Violation]:
    violations: list[Violation] = []
    legacy_public_prefixes = (
        "gelex/estimator/",
        "gelex/optim/",
        "gelex/gwas/",
        "gelex/predict/",
        "gelex/sim/",
        "gelex/logger/",
        "gelex/utils/",
        "gelex/detail/",
    )

    for file_path in collect_files(root):
        rel_file = rel_to_root(file_path, root)
        includes = iter_includes(file_path)

        for line_no, include_target in includes:
            # Rule MOD001:
            # include/gelex/internal/** only for src/** private impl
            if include_target.startswith("gelex/internal/"):
                is_internal_header = starts_with(rel_file, "include/gelex/internal")
                is_src_impl = starts_with(rel_file, "src")
                if not is_internal_header and not is_src_impl:
                    violations.append(
                        Violation(
                            rule_id="MOD001",
                            file=rel_file,
                            line=line_no,
                            include=include_target,
                            message=(
                                "`gelex/internal/**` can only be included by "
                                "`src/**` private implementation files"
                            ),
                        )
                    )

            # Rule MOD002:
            # apps/tests/benchmark must not include src/** private headers
            if (
                starts_with(rel_file, "apps")
                or starts_with(rel_file, "tests")
                or starts_with(rel_file, "benchmark")
            ) and is_rel_src_include(include_target):
                violations.append(
                    Violation(
                        rule_id="MOD002",
                        file=rel_file,
                        line=line_no,
                        include=include_target,
                        message=(
                            "`apps/tests/benchmark` must not include "
                            "`src/**` private headers"
                        ),
                    )
                )

            # Rule MOD003:
            # include/gelex/model/** must not include gelex/pipeline/**
            if starts_with(rel_file, "include/gelex/model") and include_target.startswith(
                "gelex/pipeline/"
            ):
                violations.append(
                    Violation(
                        rule_id="MOD003",
                        file=rel_file,
                        line=line_no,
                        include=include_target,
                        message="`model` public headers must not include `gelex/pipeline/**`",
                    )
                )

            # Rule MOD004:
            # io headers/impl must not include gelex/data/** or gelex/internal/data/**
            if starts_with(rel_file, "include/gelex/io") or starts_with(rel_file, "src/io"):
                if include_target.startswith("gelex/data/") or include_target.startswith(
                    "gelex/internal/data/"
                ):
                    violations.append(
                        Violation(
                            rule_id="MOD004",
                            file=rel_file,
                            line=line_no,
                            include=include_target,
                            message=(
                                "`io` headers/implementation must not include "
                                "`gelex/data/**` or `gelex/internal/data/**`"
                            ),
                        )
                    )

            # Rule MOD005:
            # types headers must not include upper-layer modules
            if starts_with(rel_file, "include/gelex/types"):
                if (
                    include_target.startswith("gelex/data/")
                    or include_target.startswith("gelex/model/")
                    or include_target.startswith("gelex/algo/")
                    or include_target.startswith("gelex/pipeline/")
                ):
                    violations.append(
                        Violation(
                            rule_id="MOD005",
                            file=rel_file,
                            line=line_no,
                            include=include_target,
                            message=(
                                "`types` headers must not include upper-layer "
                                "module headers"
                            ),
                        )
                    )

            # Rule MOD006:
            # include/gelex/** public headers must not include gelex/internal/**
            if starts_with(rel_file, "include/gelex") and not starts_with(
                rel_file, "include/gelex/internal"
            ):
                if include_target.startswith("gelex/internal/"):
                    violations.append(
                        Violation(
                            rule_id="MOD006",
                            file=rel_file,
                            line=line_no,
                            include=include_target,
                            message=(
                                "public headers under `include/gelex/**` must not "
                                "include `gelex/internal/**`"
                            ),
                        )
                    )

            # Rule MOD007:
            # legacy public include prefixes are forbidden after direct migration
            if include_target.startswith(legacy_public_prefixes):
                violations.append(
                    Violation(
                        rule_id="MOD007",
                        file=rel_file,
                        line=line_no,
                        include=include_target,
                        message=(
                            "legacy public include prefix is forbidden; use module "
                            "layout under `gelex/{infra,types,io,data,model,algo,pipeline}/**`"
                        ),
                    )
                )

            if include_target == "gelex/logger.h":
                violations.append(
                    Violation(
                        rule_id="MOD007",
                        file=rel_file,
                        line=line_no,
                        include=include_target,
                        message=(
                            "legacy public include `gelex/logger.h` is forbidden; "
                            "use `gelex/infra/logger.h`"
                        ),
                    )
                )

            if include_target in (
                "gelex/types/logger.h",
                "gelex/types/utils/math_utils.h",
                "gelex/types/utils/running_stats.h",
            ):
                replacements = {
                    "gelex/types/logger.h": "gelex/infra/logger.h",
                    "gelex/types/utils/math_utils.h": "gelex/infra/utils/math_utils.h",
                    "gelex/types/utils/running_stats.h": "gelex/infra/utils/running_stats.h",
                }
                violations.append(
                    Violation(
                        rule_id="MOD007",
                        file=rel_file,
                        line=line_no,
                        include=include_target,
                        message=(
                            "legacy public include path is forbidden; use "
                            f"`{replacements[include_target]}`"
                        ),
                    )
                )

            # Rule MOD009:
            # infra headers must not include upper-layer modules
            if starts_with(rel_file, "include/gelex/infra"):
                if (
                    include_target.startswith("gelex/types/")
                    or include_target.startswith("gelex/io/")
                    or include_target.startswith("gelex/data/")
                    or include_target.startswith("gelex/model/")
                    or include_target.startswith("gelex/algo/")
                    or include_target.startswith("gelex/pipeline/")
                ):
                    violations.append(
                        Violation(
                            rule_id="MOD009",
                            file=rel_file,
                            line=line_no,
                            include=include_target,
                            message=(
                                "`infra` headers must not include upper-layer "
                                "module headers"
                            ),
                        )
                    )

    # Rule MOD008:
    # include/gelex top-level directories must be exactly seven modules
    include_gelex = root / "include" / "gelex"
    if include_gelex.exists():
        expected_dirs = {
            "infra",
            "types",
            "io",
            "data",
            "model",
            "algo",
            "pipeline",
        }
        actual_dirs = {entry.name for entry in include_gelex.iterdir() if entry.is_dir()}
        extra_dirs = sorted(actual_dirs - expected_dirs)
        for dir_name in extra_dirs:
            violations.append(
                Violation(
                    rule_id="MOD008",
                    file="include/gelex",
                    line=1,
                    include=dir_name,
                    message=(
                        "unexpected top-level public include directory; only "
                        "`infra/types/io/data/model/algo/pipeline` are allowed"
                    ),
                )
            )

    return violations


def print_text(violations: list[Violation]) -> None:
    if not violations:
        print("No boundary violations found.")
        return

    files = {v.file for v in violations}
    print(f"Found {len(violations)} boundary violations in {len(files)} files.")
    for v in sorted(violations, key=lambda x: (x.file, x.line, x.rule_id)):
        print(
            f"- [{v.rule_id}] {v.file}:{v.line} includes `{v.include}` -> {v.message}"
        )


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]

    parser = argparse.ArgumentParser(
        description="Check module boundary include rules for Gelex."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=default_root,
        help="Repository root path (default: detected from script location).",
    )
    parser.add_argument(
        "--format",
        choices=["text", "json"],
        default="text",
        help="Output format.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()

    violations = check_boundaries(root)

    if args.format == "json":
        payload = {
            "ok": len(violations) == 0,
            "count": len(violations),
            "violations": [asdict(v) for v in violations],
        }
        print(json.dumps(payload, indent=2, ensure_ascii=False))
    else:
        print_text(violations)

    return 0 if not violations else 1


if __name__ == "__main__":
    sys.exit(main())
