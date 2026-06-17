#!/usr/bin/env python3
"""
gen_apply2Into_combinations.py
Generate the full test-combinations array for matrix.apply2Into:

  python3 gen_apply2Into_combinations.py            # Zig array to stdout
  python3 gen_apply2Into_combinations.py --stats    # also print statistics

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

Free / fixed / tied parameters
------------------------------
  * layout is always free and independent for every operand.
  * direction (permutation only) is always free and independent for every
    operand, no guard in apply2Into touches it.
  * uplo is free and independent for symmetric/hermitian operands, in any
    position, with no cross-operand constraint.
  * When O is triangular:
      - O.diag is fixed to .non_unit. This matches the
        `meta.diagOf(O) == .unit` guard in apply2Into, which always errors
        on a unit-diagonal output regardless of X/Y.
      - O.uplo is free and acts as the anchor: a triangular X or Y has its
        uplo tied to O.uplo (matches the `meta.uploOf(O) != meta.uploOf(X
        or Y)` guard). Its own diag stays free.
  * When O is not triangular, a triangular X or Y is fully free on uplo,
    diag, and layout, the uplo/diag guards in apply2Into are only ever
    emitted when O itself is triangular, so nothing is tied or fixed here.
  * When O is hermitian, any symmetric/diagonal/numeric X or Y must use a
    real element type, matching the `isComplex(...)` guards in apply2Into
    that forbid a complex symmetric/diagonal/numeric operand feeding a
    hermitian result.

N (element type) selection
----------------------------
  * X and Y always use one fixed precision, DEFAULT_INPUT_N (zsl.cf64),
    except where the hermitian-realness rule above swaps it for
    REAL_INPUT_N.
  * O normally uses DEFAULT_OUTPUT_N (zsl.cf32). Additionally, whenever
    aliasing is structurally possible, i.e., when O's variant is dense or
    static, and O's MatrixType equals X's or Y's MatrixType, a second full pass
    is emitted with O using ALIASING_OUTPUT_N (zsl.cf64). Within that
    pass, rows where the remaining free parameters also happen to
    coincide give O's pointee type and the aliased input's type genuine,
    literal equality; the other rows in that pass are simply extra
    coverage. Sparse outputs (including builder_sparse) are never
    aliasing-eligible.
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
    "builder_sparse": "matbuispa",
    "numeric": "num",
}

# Truly-dense types that disqualify sparse outputs. Diagonal and permutation
# matrices are O(n)/O(1), i.e., always sparse-compatible.
DENSE: frozenset[str] = frozenset(
    {
        "general_static",
        "general_dense",
        "symmetric_static",
        "symmetric_dense",
        "hermitian_static",
        "hermitian_dense",
        "triangular_static",
        "triangular_dense",
    }
)


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


def decide(o: str, x: str, y: str) -> str:
    """
    Returns one of:
        'unr'         unreachable
        'err:<msg>'   @compileError
        'imp:<path>'  return @import(path).apply2Into(…)
    """
    ob = base(o)

    # Invalid combinations
    if o == "numeric":
        return "unr"  # scalars cannot be matrix outputs
    if x == "builder_sparse" or y == "builder_sparse":
        return "unr"  # builder can only be output
    if ob == "permutation":
        return "err:permutation matrices cannot be used as output"
    if x == "numeric" and y == "numeric":
        return "unr"  # at least one operand must be a matrix

    # Structural output constraints (most-specific error first)
    if ob == "diagonal":
        if (x != "numeric" and base(x) != "diagonal") or (
            y != "numeric" and base(y) != "diagonal"
        ):
            return "err:diagonal output requires diagonal or numeric inputs"

    elif ob == "symmetric":
        ok = frozenset({"symmetric", "diagonal"})
        if (x != "numeric" and base(x) not in ok) or (
            y != "numeric" and base(y) not in ok
        ):
            return (
                "err:symmetric output requires symmetric, diagonal, or numeric inputs"
            )

    elif ob == "hermitian":
        # Structural check only: the real-value element-type sub-constraint is
        # handled by inline comptime guards generated by comptime_guards()
        # below.
        ok = frozenset({"hermitian", "symmetric", "diagonal"})
        if (x != "numeric" and base(x) not in ok) or (
            y != "numeric" and base(y) not in ok
        ):
            return (
                "err:hermitian output requires hermitian, symmetric, "
                "diagonal, or numeric inputs"
            )

    elif ob == "triangular":
        # Structural check only: uplo + unit-diagonal constraints are handled by
        # inline comptime guards generated by comptime_guards() below.
        ok = frozenset({"triangular", "diagonal"})
        if (x != "numeric" and base(x) not in ok) or (
            y != "numeric" and base(y) not in ok
        ):
            return (
                "err:triangular output requires triangular, diagonal, or numeric inputs"
            )

    # Sparse-output density constraint
    if variant(o) == "sparse" and (x in DENSE or y in DENSE):
        return "err:sparse output requires sparse or numeric inputs (no dense operand)"

    return f"imp:apply2/{FID[o]}_{FID[x]}_{FID[y]}.zig"


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

DEFAULT_OUTPUT_N = "zsl.cf32"
ALIASING_OUTPUT_N = "zsl.cf64"
DEFAULT_INPUT_N = "zsl.cf64"
REAL_INPUT_N = "f64"

Role = str  # "O" | "X" | "Y"
Axis = tuple[Role, str]  # e.g. ("X", "uplo")


def operand_axis_info(
    role: Role, tv: str, o_kind: str
) -> tuple[list[tuple[str, list[str]]], dict[str, str]]:
    """Free axes (name, candidate values) and fixed values for one operand,
    not counting N (which is handled separately by operand_n)."""
    if tv == "numeric":
        return [], {}
    kind = base(tv)
    if kind == "general":
        return [("layout", LAYOUTS)], {}
    if kind in ("symmetric", "hermitian"):
        return [("uplo", UPLOS), ("layout", LAYOUTS)], {}
    if kind in ("diagonal", "builder_sparse"):
        return [], {}
    if kind == "permutation":
        return [("direction", DIRECTIONS)], {}
    if kind == "triangular":
        if role == "O":
            return [("uplo", UPLOS), ("layout", LAYOUTS)], {"diag": "non_unit"}
        if o_kind == "triangular":
            return [("diag", DIAGS), ("layout", LAYOUTS)], {}  # uplo tied to O's anchor
        return [("uplo", UPLOS), ("diag", DIAGS), ("layout", LAYOUTS)], {}
    raise AssertionError(f"unhandled kind {kind!r} for type {tv!r}")


def needs_real(o_kind: str, tv: str) -> bool:
    """True if this X/Y operand must use a real (non-complex) element type
    because O is hermitian and this operand is symmetric, diagonal, or
    bare numeric."""
    if o_kind != "hermitian":
        return False
    if tv == "numeric":
        return True
    return base(tv) in ("symmetric", "diagonal")


def input_n(tv: str, o_kind: str) -> str:
    return REAL_INPUT_N if needs_real(o_kind, tv) else DEFAULT_INPUT_N


def render_operand(
    tv: str, n_expr: str, role: Role, assignment: dict[Axis, str]
) -> str:
    if tv == "numeric":
        return n_expr
    kind = base(tv)
    var = variant(tv)
    prefix = "tzsl" if var == "static" else "zsl"
    ctor = CTOR[var]
    ns = KIND_NS[kind]
    params = [n_expr]
    if kind == "permutation":
        params.append(f".{assignment[(role, 'direction')]}")
    if kind in ("symmetric", "hermitian", "triangular"):
        params.append(f".{assignment[(role, 'uplo')]}")
    if kind == "triangular":
        params.append(f".{assignment[(role, 'diag')]}")
    if kind in ("general", "symmetric", "hermitian", "triangular"):
        params.append(f".{assignment[(role, 'layout')]}")
    return f"{prefix}.matrix.{ns}.{ctor}(" + ", ".join(params) + ")"


def generate_group(o: str, x: str, y: str) -> tuple[str, list[str]]:
    """Comment + every concrete `.{ O, X, Y },` row for one valid MatrixType
    combination (all free-axis values crossed, doubled for aliasing where
    eligible)."""
    o_kind = base(o)

    axis_specs: list[tuple[Role, str, list[str]]] = []
    fixed: dict[Axis, str] = {}
    tied: dict[Axis, Axis] = {}

    for role, tv in (("O", o), ("X", x), ("Y", y)):
        free, fx = operand_axis_info(role, tv, o_kind)
        for axis_name, values in free:
            axis_specs.append((role, axis_name, values))
        for axis_name, val in fx.items():
            fixed[(role, axis_name)] = val
        if (
            role in ("X", "Y")
            and tv != "numeric"
            and base(tv) == "triangular"
            and o_kind == "triangular"
        ):
            tied[(role, "uplo")] = ("O", "uplo")

    names = [(r, a) for (r, a, _) in axis_specs]
    value_lists = [vals for (_, _, vals) in axis_specs]

    x_n = input_n(x, o_kind)
    y_n = input_n(y, o_kind)

    aliasing_eligible = variant(o) in ("dense", "static") and (o == x or o == y)
    o_n_list = (
        [DEFAULT_OUTPUT_N, ALIASING_OUTPUT_N]
        if aliasing_eligible
        else [DEFAULT_OUTPUT_N]
    )

    lines: list[str] = []
    for o_n in o_n_list:
        for combo_values in itertools.product(*value_lists):
            assignment: dict[Axis, str] = dict(zip(names, combo_values))
            assignment.update(fixed)
            for k, anchor in tied.items():
                assignment[k] = assignment[anchor]

            o_expr = render_operand(o, o_n, "O", assignment)
            x_expr = render_operand(x, x_n, "X", assignment)
            y_expr = render_operand(y, y_n, "Y", assignment)
            lines.append(f".{{ {o_expr}, {x_expr}, {y_expr} }},")

    return f"// {FID[o]}_{FID[x]}_{FID[y]}", lines


def all_groups() -> list[tuple[str, str, str, str, list[str]]]:
    """Every valid (o, x, y) triple with its rendered comment + rows."""
    out = []
    for o in TYPES:
        for x in TYPES:
            for y in TYPES:
                if not decide(o, x, y).startswith("imp:"):
                    continue
                comment, lines = generate_group(o, x, y)
                out.append((o, x, y, comment, lines))
    return out


def generate(
    base_indent: int = 8,
) -> tuple[str, list[tuple[str, str, str, str, list[str]]]]:
    bi = " " * base_indent
    groups = all_groups()
    total_lines = sum(len(lines) for *_, lines in groups)
    quota = max(50_000, total_lines * 50 + 100_000)

    out = [
        "const combinations = blk: {",
        f"    @setEvalBranchQuota({quota});",
        "    break :blk [_][3]type{",
    ]
    for i, (_, _, _, comment, lines) in enumerate(groups):
        if i > 0:
            out.append("")
        out.append(f"{bi}{comment}")
        out.extend(f"{bi}{ln}" for ln in lines)
    out.append("    };")
    out.append("};")
    return "\n".join(out), groups


if __name__ == "__main__":
    text, groups = generate()

    if "--stats" in sys.argv:
        total_lines = sum(len(lines) for *_, lines in groups)
        doubled = sum(
            1
            for o, x, y, _, _ in groups
            if variant(o) in ("dense", "static") and (o == x or o == y)
        )
        real_needed = sum(
            1
            for o, x, y, _, _ in groups
            if needs_real(base(o), x) or needs_real(base(o), y)
        )
        sizes = Counter(len(lines) for *_, lines in groups)
        print(f"Valid groups          : {len(groups):>7}", file=sys.stderr)
        print(f"Total rows            : {total_lines:>7}", file=sys.stderr)
        print(f"Aliasing-doubled grps : {doubled:>7}", file=sys.stderr)
        print(f"Real-input-N groups   : {real_needed:>7}", file=sys.stderr)
        print(
            f"Smallest/largest group: {min(sizes)} / {max(sizes)} rows", file=sys.stderr
        )

    print(text)
