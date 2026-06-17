#!/usr/bin/env python3
"""
gen_apply2_combinations.py
Generate the full test-combinations array for matrix.apply2Alloc:

  python3 gen_apply2_combinations.py            # Zig array to stdout
  python3 gen_apply2_combinations.py --stats    # also print statistics

Type construction
------------------
  zsl.matrix.<kind>.<Static|Dense|Sparse>(...)
    general                -> N, layout
    symmetric, hermitian   -> N, uplo, layout
    triangular             -> N, uplo, diag, layout
    permutation            -> N, direction
    diagonal, builder      -> N
  numeric operand          -> bare N
A *_static variant uses the "tzsl" prefix instead of "zsl" and omits the
dimension parameters); every other parameter keeps its normal position.
"""

from __future__ import annotations

import itertools
import sys
from collections import Counter

# MatrixType enum members
TYPES: list[str] = [
    "general_static",
    "general_dense",
    "general_sparse",
    "symmetric_static",
    "symmetric_dense",
    "symmetric_sparse",
    "hermitian_static",
    "hermitian_dense",
    "hermitian_sparse",
    "triangular_static",
    "triangular_dense",
    "triangular_sparse",
    "diagonal_static",
    "diagonal_sparse",
    "permutation_static",
    "permutation_sparse",
    "builder_sparse",
    "numeric",
]

# File-name id fragments (mat<type><variant>  or  num):
#   general    -> gen    symmetric  -> sym    hermitian   -> her
#   triangular -> tri    diagonal   -> dia    permutation -> per    builder -> bui
#   static     -> sta    dense      -> den    sparse      -> spa
FID: dict[str, str] = {
    "general_static": "matgensta",
    "general_dense": "matgenden",
    "general_sparse": "matgenspa",
    "symmetric_static": "matsymsta",
    "symmetric_dense": "matsymden",
    "symmetric_sparse": "matsymspa",
    "hermitian_static": "mathersta",
    "hermitian_dense": "matherden",
    "hermitian_sparse": "matherspa",
    "triangular_static": "mattrista",
    "triangular_dense": "mattriden",
    "triangular_sparse": "mattrispa",
    "diagonal_static": "matdiasta",
    "diagonal_sparse": "matdiaspa",
    "permutation_static": "matpersta",
    "permutation_sparse": "matperspa",
    # "builder_sparse": "matbuispa",
    "numeric": "num",
}


def base(t: str) -> str:
    """'hermitian_dense' -> 'hermitian';  'numeric'/'builder_sparse' -> itself."""
    if t in ("numeric", "builder_sparse"):
        return t
    return "_".join(t.split("_")[:-1])


def variant(t: str) -> str | None:
    """'triangular_sparse' -> 'sparse'; 'numeric' -> None."""
    if t == "numeric":
        return None
    return t.rsplit("_", 1)[-1]


LAYOUTS = ["row_major", "col_major"]
UPLOS = ["upper", "lower"]
DIAGS = ["non_unit", "unit"]
DIRECTIONS = ["forward", "backward"]

KIND_NS = {
    "general": "general",
    "symmetric": "symmetric",
    "hermitian": "hermitian",
    "triangular": "triangular",
    "diagonal": "diagonal",
    "permutation": "permutation",
    "builder_sparse": "builder",
}
CTOR = {"static": "Static", "dense": "Dense", "sparse": "Sparse"}

INPUT_N = "zsl.cf64"


def operand_axes(tv: str) -> list[tuple[str, list[str]]]:
    """Free (axis_name, candidate values) for one operand, not counting N."""
    if tv == "numeric":
        return []
    kind = base(tv)
    if kind == "general":
        return [("layout", LAYOUTS)]
    if kind in ("symmetric", "hermitian"):
        return [("uplo", UPLOS), ("layout", LAYOUTS)]
    if kind == "triangular":
        return [("uplo", UPLOS), ("diag", DIAGS), ("layout", LAYOUTS)]
    if kind == "permutation":
        return [("direction", DIRECTIONS)]
    if kind in ("diagonal", "builder_sparse"):
        return []
    raise AssertionError(f"unhandled kind {kind!r} for type {tv!r}")


def render_operand(tv: str, values: dict[str, str]) -> str:
    if tv == "numeric":
        return INPUT_N
    kind = base(tv)
    var = variant(tv)
    prefix = "tzsl" if var == "static" else "zsl"
    ctor = CTOR[var]
    ns = KIND_NS[kind]
    params = [INPUT_N]
    if kind == "permutation":
        params.append(f".{values['direction']}")
    if kind in ("symmetric", "hermitian", "triangular"):
        params.append(f".{values['uplo']}")
    if kind == "triangular":
        params.append(f".{values['diag']}")
    if kind in ("general", "symmetric", "hermitian", "triangular"):
        params.append(f".{values['layout']}")
    return f"{prefix}.matrix.{ns}.{ctor}(" + ", ".join(params) + ")"


def generate_pair(x: str, y: str) -> tuple[str, list[str]]:
    """Comment + every concrete `.{ X, Y },` row for one (X, Y) MatrixType
    pair."""
    x_axes = operand_axes(x)
    y_axes = operand_axes(y)

    x_names = [a for a, _ in x_axes]
    y_names = [a for a, _ in y_axes]
    value_lists = [v for _, v in x_axes] + [v for _, v in y_axes]
    split = len(x_axes)

    lines: list[str] = []
    for combo in itertools.product(*value_lists):
        x_values = dict(zip(x_names, combo[:split]))
        y_values = dict(zip(y_names, combo[split:]))
        x_expr = render_operand(x, x_values)
        y_expr = render_operand(y, y_values)
        lines.append(f".{{ {x_expr}, {y_expr} }},")

    return f"// {FID[x]}_{FID[y]}", lines


def all_pairs() -> list[tuple[str, str, str, list[str]]]:
    """Every (x, y) MatrixType pair -- the full 18*18 grid, unfiltered."""
    return [(x, y, *generate_pair(x, y)) for x in TYPES for y in TYPES]


def generate(base_indent: int = 8) -> tuple[str, list[tuple[str, str, str, list[str]]]]:
    bi = " " * base_indent
    pairs = all_pairs()
    total_lines = sum(len(lines) for *_, lines in pairs)
    quota = max(50_000, total_lines * 50 + 100_000)

    out = [
        "const combinations = blk: {",
        f"    @setEvalBranchQuota({quota}); // heuristic -- bump if Zig complains",
        "    break :blk [_][2]type{",
    ]
    for i, (_, _, comment, lines) in enumerate(pairs):
        if i > 0:
            out.append("")
        out.append(f"{bi}{comment}")
        out.extend(f"{bi}{ln}" for ln in lines)
    out.append("    };")
    out.append("};")
    return "\n".join(out), pairs


if __name__ == "__main__":
    text, pairs = generate()

    if "--stats" in sys.argv:
        total_lines = sum(len(lines) for *_, lines in pairs)
        sizes = Counter(len(lines) for *_, lines in pairs)
        print(
            f"Kind-level pairs : {len(pairs):>7}  (expect 18*18 = 324)", file=sys.stderr
        )
        print(f"Total rows       : {total_lines:>7}", file=sys.stderr)
        print(
            f"Smallest/largest pair: {min(sizes)} / {max(sizes)} rows", file=sys.stderr
        )

    print(text)
