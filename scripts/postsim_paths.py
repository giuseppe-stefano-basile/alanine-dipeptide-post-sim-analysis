#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable


GLOB_META = set("*?[]")


class InputResolutionError(RuntimeError):
    """Raised when a configured input cannot be resolved unambiguously."""


def has_glob(pattern: str) -> bool:
    return any(ch in pattern for ch in GLOB_META)


def load_cases(cases_json: Path) -> list[dict]:
    obj = json.loads(cases_json.read_text())
    cases = obj.get("cases", [])
    if not isinstance(cases, list):
        raise InputResolutionError(f"'cases' must be a list in {cases_json}")
    return cases


def normalize_search_roots(search_roots: Iterable[str]) -> list[Path]:
    roots = []
    for raw in search_roots:
        if raw is None:
            continue
        for item in str(raw).split(":"):
            item = item.strip()
            if not item:
                continue
            root = Path(item).expanduser().resolve()
            if not root.is_dir():
                raise InputResolutionError(f"Search root is not a directory: {root}")
            roots.append(root)

    deduped = []
    seen = set()
    for root in roots:
        key = str(root)
        if key not in seen:
            deduped.append(root)
            seen.add(key)

    if not deduped:
        raise InputResolutionError(
            "No search root was provided. Choose at least one directory with --search-root."
        )
    return deduped


def patterns_from_value(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, list):
        out = []
        for item in value:
            out.extend(patterns_from_value(item))
        return out
    if isinstance(value, dict):
        out = []
        for key in (
            "path",
            "paths",
            "pattern",
            "patterns",
            "filename",
            "filenames",
            "glob",
            "globs",
        ):
            if key in value:
                out.extend(patterns_from_value(value[key]))
        return out
    raise InputResolutionError(f"Unsupported input spec: {value!r}")


def _dedupe(paths: Iterable[Path]) -> list[Path]:
    out = []
    seen = set()
    for path in paths:
        absolute = path.expanduser().absolute()
        key = str(absolute)
        if key not in seen:
            out.append(absolute)
            seen.add(key)
    return out


def _direct_matches(pattern: str, roots: list[Path]) -> list[Path]:
    path = Path(pattern).expanduser()
    if path.is_absolute():
        return [path] if path.is_file() else []

    matches = []
    for root in roots:
        if has_glob(pattern):
            matches.extend(p for p in root.glob(pattern) if p.is_file())
        else:
            candidate = root / pattern
            if candidate.is_file():
                matches.append(candidate)
    return _dedupe(matches)


def _recursive_matches(pattern: str, roots: list[Path]) -> list[Path]:
    path = Path(pattern)
    matches = []

    for root in roots:
        if has_glob(pattern):
            matches.extend(p for p in root.rglob(pattern) if p.is_file())
            continue

        parts = path.parts
        if not parts:
            continue
        basename = parts[-1]
        for candidate in root.rglob(basename):
            if not candidate.is_file():
                continue
            if len(parts) == 1:
                matches.append(candidate)
                continue
            try:
                rel_parts = candidate.relative_to(root).parts
            except ValueError:
                continue
            if tuple(rel_parts[-len(parts) :]) == tuple(parts):
                matches.append(candidate)

    return _dedupe(matches)


def search_pattern(pattern: str, roots: list[Path]) -> list[Path]:
    pattern = str(pattern).strip()
    if not pattern:
        return []

    direct = _direct_matches(pattern, roots)
    if direct:
        return direct
    return _recursive_matches(pattern, roots)


def resolve_input(
    *,
    case_name: str,
    field: str,
    value,
    roots: list[Path],
    required: bool,
) -> tuple[str, list[str]]:
    patterns = patterns_from_value(value)
    warnings = []

    if not patterns:
        if required:
            raise InputResolutionError(f"Case '{case_name}' is missing required field '{field}'")
        return "", warnings

    searched = []
    for pattern in patterns:
        searched.append(str(pattern))
        candidates = search_pattern(pattern, roots)
        if not candidates:
            continue
        if len(candidates) == 1:
            return str(candidates[0]), warnings
        shown = "\n".join(f"  - {p}" for p in candidates[:20])
        more = "" if len(candidates) <= 20 else f"\n  ... and {len(candidates) - 20} more"
        raise InputResolutionError(
            f"Case '{case_name}' field '{field}' is ambiguous for pattern '{pattern}'.\n"
            f"Choose a narrower --search-root or make the config pattern more specific:\n"
            f"{shown}{more}"
        )

    roots_txt = ", ".join(str(root) for root in roots)
    msg = (
        f"Case '{case_name}' field '{field}' was not found.\n"
        f"Patterns searched: {searched}\n"
        f"Search roots: {roots_txt}"
    )
    if required:
        raise InputResolutionError(msg)
    warnings.append(msg)
    return "", warnings


def resolve_case(case: dict, roots: list[Path], strict_logs: bool = False) -> tuple[dict, list[str]]:
    name = str(case.get("name", "")).strip()
    if not name:
        raise InputResolutionError("Every case must define a non-empty 'name'")

    warnings: list[str] = []
    resolved = {
        "name": name,
        "description": str(case.get("description", "")),
    }

    for field in ("npbc_prod_dump", "pbc_prod_dump"):
        resolved[field], field_warnings = resolve_input(
            case_name=name,
            field=field,
            value=case.get(field),
            roots=roots,
            required=True,
        )
        warnings.extend(field_warnings)

    for field in ("npbc_prod_log", "pbc_prod_log"):
        resolved[field], field_warnings = resolve_input(
            case_name=name,
            field=field,
            value=case.get(field),
            roots=roots,
            required=bool(strict_logs),
        )
        warnings.extend(field_warnings)

    return resolved, warnings


def select_cases(cases: list[dict], requested: Iterable[str] | None) -> list[dict]:
    requested_names = [name for name in (requested or []) if name]
    if not requested_names:
        return cases

    by_name = {str(case.get("name", "")): case for case in cases}
    selected = []
    missing = []
    for name in requested_names:
        if name in by_name:
            selected.append(by_name[name])
        else:
            missing.append(name)
    if missing:
        known = ", ".join(sorted(by_name))
        raise InputResolutionError(f"Unknown case(s): {', '.join(missing)}. Known cases: {known}")
    return selected
