#!/usr/bin/env python3

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import Mapping, MutableMapping

BASE_DIR = Path(__file__).resolve().parents[1]
MATRICES_DIR = BASE_DIR / "Matrices"
KNOWN_PERMANENTS_FILE = BASE_DIR / "known_permanents.md"

SECTION_ORDER = [
    "TinyOriginal",
    "Banded",
    "BlockDiagonal",
    "FromTinyToLarge",
]

SECTION_NOTES = {
    "TinyOriginal": (
        "Entries in this section come from "
        "`/home/kaya/Rasmussen/CleanBackup/PermanentResults/small_permanent`."
    ),
    "Banded": "These permanents are exact and come from the cyclic band generator.",
    "BlockDiagonal": (
        "These permanents are exact because the matrices are block diagonal and "
        "each block permanent is computed exactly."
    ),
    "FromTinyToLarge": (
        "These permanents are products of the TinyOriginal values listed in this "
        "file."
    ),
}

INTRO_LINES = [
    "# Known Permanents Database",
    "",
    "This file tracks matrices in `Matrices/` whose permanent is known or is",
    "treated as ground truth for the current experiments.",
    "",
    "Matrices with permanent `0` are intentionally omitted.",
    "Families such as `Small`, `Medium`, `Large`, and `ErdosRenyi` are not listed",
    "here because their permanents are currently unknown.",
    "",
]


def repo_relative(path: str | Path) -> str:
    path_obj = Path(path)
    if not path_obj.is_absolute():
        return path_obj.as_posix()
    return path_obj.resolve().relative_to(BASE_DIR).as_posix()


def read_matrix_market_metadata(path: str | Path) -> OrderedDict[str, str]:
    matrix_path = Path(path)
    market_header = None

    with matrix_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if market_header is None:
                if not line.startswith("%%MatrixMarket"):
                    raise ValueError(f"{matrix_path} is not a Matrix Market file")
                market_header = line.removeprefix("%%MatrixMarket ").strip()
                continue
            if line.startswith("%"):
                continue
            rows, cols, stored_entries = map(int, line.split()[:3])
            return OrderedDict(
                [
                    ("Rows", str(rows)),
                    ("Cols", str(cols)),
                    ("Stored entries", str(stored_entries)),
                    ("MatrixMarket", market_header),
                ]
            )

    raise ValueError(f"Could not read Matrix Market metadata from {matrix_path}")


def parse_registry(md_path: str | Path = KNOWN_PERMANENTS_FILE) -> OrderedDict[str, OrderedDict[str, OrderedDict[str, str]]]:
    path = Path(md_path)
    sections: OrderedDict[str, OrderedDict[str, OrderedDict[str, str]]] = OrderedDict(
        (name, OrderedDict()) for name in SECTION_ORDER
    )
    if not path.exists():
        return sections

    current_section = None
    current_entry_path = None
    current_fields: OrderedDict[str, str] | None = None

    def finalize_entry() -> None:
        nonlocal current_entry_path, current_fields
        if current_section and current_entry_path and current_fields is not None:
            sections.setdefault(current_section, OrderedDict())[current_entry_path] = current_fields
        current_entry_path = None
        current_fields = None

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()

            if stripped.startswith("## "):
                finalize_entry()
                current_section = stripped[3:].strip()
                sections.setdefault(current_section, OrderedDict())
                continue

            if stripped.startswith("### "):
                finalize_entry()
                current_entry_path = stripped[4:].strip()
                current_fields = OrderedDict()
                continue

            if not stripped or stripped == "No entries yet.":
                continue

            if current_entry_path and ":" in stripped:
                key, value = stripped.split(":", 1)
                current_fields[key.strip()] = value.strip()

    finalize_entry()
    return sections


def write_registry(
    sections: MutableMapping[str, MutableMapping[str, Mapping[str, str]]],
    md_path: str | Path = KNOWN_PERMANENTS_FILE,
) -> None:
    ordered_sections = list(SECTION_ORDER)
    for section_name in sections:
        if section_name not in ordered_sections:
            ordered_sections.append(section_name)

    lines = []
    for line in INTRO_LINES:
        if line:
            lines.append(f"{line}\n")
        else:
            lines.append("\n")
    for section_name in ordered_sections:
        lines.append(f"## {section_name}\n")
        note = SECTION_NOTES.get(section_name)
        if note:
            lines.append(f"{note}\n")
        lines.append("\n")

        entries = sections.get(section_name, OrderedDict())
        if not entries:
            lines.append("No entries yet.\n\n")
            continue

        for relative_path in sorted(entries):
            lines.append(f"### {relative_path}\n")
            for key, value in entries[relative_path].items():
                if value is None or value == "":
                    continue
                lines.append(f"{key}: {value}\n")
            lines.append("\n")

    Path(md_path).write_text("".join(lines))


def upsert_registry_entry(
    section: str,
    relative_path: str | Path,
    fields: Mapping[str, object],
    md_path: str | Path = KNOWN_PERMANENTS_FILE,
) -> None:
    sections = parse_registry(md_path)
    section_entries = sections.setdefault(section, OrderedDict())
    normalized_fields = OrderedDict()
    for key, value in fields.items():
        if value is None:
            continue
        text = str(value)
        if text == "":
            continue
        normalized_fields[key] = text
    section_entries[repo_relative(relative_path)] = normalized_fields
    write_registry(sections, md_path)
