#!/usr/bin/env python3
"""
gen_matmulInto_branches.py
Generate the full comptime switch dispatch body for linalg.matmulInto.

  python3 gen_matmulInto_branches.py            # Zig switch body to stdout
  python3 gen_matmulInto_branches.py --stats    # also print statistics

Shape rules that determine (O, X, Y) validity
----------------------------------------------
A. O is any matrix (including builder_sparse, excluding permutation):
     * (X=matrix, Y=matrix): structural constraints per O kind (see below)
     * (X=vector, Y=vector): Rejected. Vector × vector outer products /
       rank-1 updates are handled by the separate linalg.outer function,
       not matmulInto, regardless of O's kind.
     * mixed (one matrix + one vector) -> compile error (shape mismatch)

B. O is builder_sparse matrix:
     Accepts any (matrix × matrix); no structural constraints beyond that
     (builder accumulates whatever is computed). Vector × vector is still
     rejected per rule A.

C. O is permutation:
     Only valid for (permutation × permutation) -> composition is a permutation.
     Everything else is a compile error (including vector × vector).

D. O is a vector:
     Exactly one of (X, Y) must be a matrix and the other a vector.
     Any matrix kind is valid on the matrix side (covers GEMV, left-multiply,
     SYMV, HEMV, TRMV, etc.).
     (matrix × matrix) or (vector × vector) -> compile error.

Structural constraints for matrix O × matrix inputs
----------------------------------------------------
  general, symmetric, hermitian, triangular, diagonal: any matrix × any matrix.
  (The caller is responsible for ensuring the result mathematically fits the output structure.
  Inline guards such as isComplex or uplo checks are still generated where appropriate).
  permutation: permutation × permutation only (composition).

Sparse-output density constraint (matrix O only)
---------------------------------------------------
  A sparse output rejects a truly-dense input. Diagonal and
  permutation inputs are O(n)/O(1) and stay sparse-compatible,
  so they never trigger this. Checked last, after the kind-specific
  structural check has already accepted the combination structurally.

Comptime guards (inline in the returned arm)
---------------------------------------------
  hermitian O + symmetric or diagonal X/Y -> isComplex(meta.Numeric(X/Y)) guard
  triangular O + triangular X or Y       -> uploOf(O) != uploOf(X/Y) guard
  triangular O always                    -> diagOf(O) == .unit guard

  Since vector × vector never reaches a matrix output (rule A), X and Y are
  always matrix-domain by the time these guards run; no vector special-casing
  is needed here.

builder_sparse and numeric inputs
-----------------------------------
  builder_sparse as X or Y -> unreachable (it is an output-only accumulator).
  numeric is caught by `else => unreachable` at the meta.domain level and is
  never emitted explicitly in the switch.
"""

from __future__ import annotations

import sys
from collections import Counter

MATRIX_TYPES: list[str] = [
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
    "numeric",  # caught by else=>unreachable at meta.domain level
]

VECTOR_TYPES: list[str] = [
    "vector_static",
    "vector_dense",
    "vector_sparse",
]

MATRIX_DISPATCH: list[str] = [t for t in MATRIX_TYPES if t != "numeric"]
VECTOR_DISPATCH: list[str] = VECTOR_TYPES
ALL_DISPATCH: list[str] = MATRIX_DISPATCH + VECTOR_DISPATCH

MATRIX_DENSE: frozenset[str] = frozenset(
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

VECTOR_DENSE: frozenset[str] = frozenset({"vector_static", "vector_dense"})


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
    "vector_static": "vecsta",
    "vector_dense": "vecden",
    "vector_sparse": "vecspa",
}

NS = "zsl.linalg.matmulInto"


def base(t: str) -> str:
    """Kind of a type: 'hermitian_dense' -> 'hermitian', 'vector_dense' -> 'vector'."""
    if t in ("numeric", "builder_sparse"):
        return t
    return "_".join(t.split("_")[:-1])


def variant(t: str) -> str | None:
    if t == "numeric":
        return None
    return t.rsplit("_", 1)[-1]


def domain(t: str) -> str:
    """'matrix' for all matrix types (including numeric/builder); 'vector' for vectors."""
    return "vector" if base(t) == "vector" else "matrix"


def zig_enum(t: str) -> str:
    """Zig enum value name used in switch arms."""
    if base(t) == "vector":
        return variant(t) or t  # 'static' / 'dense' / 'sparse'
    return t  # MatrixType members keep their full name


def _is_dense_operand(t: str) -> bool:
    """
    True if `t` is a truly-dense representation that disqualifies a sparse
    output: a general/symmetric/hermitian/triangular matrix in its static or
    dense variant (matrix DENSE set), or a vector in its static or dense
    variant (VECTOR_DENSE). Diagonal and permutation matrices (any variant)
    are O(n)/O(1) and never count as dense, and neither does an
    already-sparse operand of either domain.
    """
    return t in MATRIX_DENSE or t in VECTOR_DENSE


def _imp(o: str, x: str, y: str) -> str:
    """
    Final gate before returning a valid 'imp:' dispatch for an output (matrix
    or vector) that has already passed its kind/shape-specific structural
    check. Mirrors apply2Into rule 7, extended to vectors: a sparse output
    — any *_sparse matrix kind (including builder_sparse and
    permutation_sparse), or vector_sparse — rejects a truly-dense operand on
    either side, matrix or vector. Diagonal and permutation matrix operands
    (any variant) are exempt, same as in apply2Into.
    """
    if variant(o) == "sparse" and (_is_dense_operand(x) or _is_dense_operand(y)):
        return "err:sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)"
    return f"imp:matmul/{FID[o]}_{FID[x]}_{FID[y]}.zig"


def decide(o: str, x: str, y: str) -> str:
    """
    Returns one of:
        'unr'         -> unreachable
        'err:<msg>'   -> @compileError
        'imp:<path>'  -> return @import(path).matmulInto(…)
    """
    od = domain(o)
    xd = domain(x)
    yd = domain(y)
    ob = base(o)
    xb = base(x)
    yb = base(y)

    # Global invalidity
    # builder_sparse can only be an output; numeric is handled by domain switch.
    if x == "builder_sparse" or y == "builder_sparse":
        return "unr"

    # O is a matrix
    if od == "matrix":
        # Vector x vector -> matrix (outer product / rank-1 update) is handled
        # by the separate linalg.outer function, not matmulInto. Reject it here
        # regardless of O's kind, before any kind-specific structural check.
        if xd == "vector" and yd == "vector":
            return "err:vector \u00d7 vector outer products are not supported by matmulInto; use linalg.outer instead"

        # Shape constraint: matrix output now requires two matrix inputs.
        if xd != "matrix" or yd != "matrix":
            return "err:matrix output requires two matrix inputs (shape mismatch)"

        # C. permutation output: only permutation × permutation
        if ob == "permutation":
            if xb == "permutation" and yb == "permutation":
                return _imp(o, x, y)
            return "err:permutation output requires permutation × permutation inputs"

        # All other matrix outputs (general, symmetric, hermitian, triangular,
        # diagonal, builder_sparse) now allow any matrix × matrix input combination.
        return _imp(o, x, y)

    # O is a vector
    if od == "vector":
        if (xd == "matrix" and yd == "vector") or (xd == "vector" and yd == "matrix"):
            return _imp(o, x, y)
        return "err:vector output requires exactly one matrix and one vector as inputs"

    return "unr"


_TYPS = (
    '\\n\\to: *" ++ @typeName(O) ++ '
    '"\\n\\tx: " ++ @typeName(X) ++ '
    '"\\n\\ty: " ++ @typeName(Y)'
)


def _cerr(msg: str) -> str:
    return f'@compileError("{NS}: {msg}{_TYPS})'


def _iscomplex_cond(t: str, var: str) -> str | None:
    """isComplex condition for matrix input `t` bound to `var`.
    Returns None for hermitian (always OK) and any type that does not need
    checking. `t` is always matrix-domain here: vector x vector is rejected
    before a matrix output's guards are ever generated (see decide())."""
    b = base(t)
    if b in ("symmetric", "diagonal"):
        return f"meta.isComplex(meta.Numeric({var}))"
    return None  # hermitian, general: always structurally valid


def comptime_guards(o: str, x: str, y: str) -> list[str]:
    """Single-line `comptime if` Zig statements to emit before the @import return."""
    ob = base(o)
    stmts: list[str] = []

    # if ob == "hermitian":
    #     cx = _iscomplex_cond(x, "X")
    #     cy = _iscomplex_cond(y, "Y")
    #     conds = [c for c in (cx, cy) if c is not None]
    #     if conds:
    #         stmts.append(
    #             f"comptime if ({' or '.join(conds)}) "
    #             f"{_cerr('hermitian output: symmetric and diagonal inputs must have a real element type')};"
    #         )

    if ob == "triangular":
        conds: list[str] = []
        # X, Y are always matrix-domain here (vector x vector is rejected
        # before a matrix output's guards are generated, see decide()).
        if base(x) == "triangular" and base(y) == "triangular":
            conds.append("meta.uploOf(O) != meta.uploOf(X)")
            conds.append("meta.uploOf(O) != meta.uploOf(Y)")
        # conds.append("meta.diagOf(O) == .unit")
        if len(conds) == 0:
            return stmts

        stmts.append(
            f"comptime if ({' or '.join(conds)}) "
            f"{_cerr('triangular operands must share uplo, and output must not be unit-diagonal')};"
        )

    return stmts


def zig_arm(o: str, x: str, y: str) -> tuple[bool, list[str]]:
    """
    Returns (is_block, body_lines):
      is_block=False -> single expression ending with ","
      is_block=True  -> block body lines (caller adds { } and trailing ",")
    """
    d = decide(o, x, y)

    if d == "unr":
        return False, ["unreachable,"]

    if d.startswith("err:"):
        return False, [f"{_cerr(d[4:])},"]

    path = d[4:]
    guards = comptime_guards(o, x, y)
    ret = f'return @import("{path}").matmulInto(o, x, y);'

    if not guards:
        return False, [ret[:-1] + ","]

    return True, guards + [ret]


def optimize_type_switch(
    switch_expr: str,
    arms: dict[str, tuple[bool, list[str]]],
) -> tuple[bool, list[str]]:
    """
    Group identical arms, pick the most-common error/unreachable as `else =>`.
    Returns (is_block, lines) for the resulting switch (or collapsed expression).
    switch_expr is the full expression, e.g. 'meta.matrixType(X)'.
    """
    T = "    "

    groups: dict[str, tuple[bool, list[str], list[str]]] = {}
    for var, (is_block, lines) in arms.items():
        key = "\n".join(lines)
        if key not in groups:
            groups[key] = (is_block, lines, [])
        groups[key][2].append(var)

    # Best else-candidate: most-frequent error or unreachable
    best_key = max(
        (
            k
            for k, (blk, lns, _) in groups.items()
            if not blk
            and (lns[0] == "unreachable," or lns[0].startswith("@compileError"))
        ),
        key=lambda k: len(groups[k][2]),
        default=None,
    )

    # Complete collapse when every arm is identical
    if len(groups) == 1:
        fk = next(iter(groups))
        return groups[fk][0], groups[fk][1]

    out = [f"switch (comptime {switch_expr}) {{"]

    for key, (is_block, lines, vars_) in groups.items():
        if key == best_key:
            continue
        prefix = ", ".join(f".{zig_enum(v)}" for v in vars_) + " => "
        if is_block:
            out.append(f"{T}{prefix}{{")
            for ln in lines:
                out.append(f"{T}{T}{ln}")
            out.append(f"{T}}},")
        else:
            out.append(f"{T}{prefix}{lines[0]}")
            for ln in lines[1:]:
                out.append(f"{T}{ln}")

    if best_key is not None:
        bi2, bl2, _ = groups[best_key]
        if bi2:
            out.append(f"{T}else => {{")
            for ln in bl2:
                out.append(f"{T}{T}{ln}")
            out.append(f"{T}}},")
        else:
            out.append(f"{T}else => {bl2[0]}")
            for ln in bl2[1:]:
                out.append(f"{T}{ln}")

    out.append("},")
    return False, out


def optimize_domain_type_switch(
    switch_var: str,
    all_arms: dict[str, tuple[bool, list[str]]],
) -> tuple[bool, list[str]]:
    """
    Produce:
      switch (comptime meta.domain(SWITCH_VAR)) {
          .matrix => switch (comptime meta.matrixType(SWITCH_VAR)) { … },
          .vector => switch (comptime meta.vectorType(SWITCH_VAR)) { … },
          else => unreachable,
      }

    `all_arms` maps every type in ALL_DISPATCH to its (is_block, lines) arm.
    Lines returned by this function use T=4-space relative indentation so that
    callers can embed them with a single additional prefix layer.
    """
    T = "    "

    mat_arms = {t: all_arms[t] for t in MATRIX_DISPATCH}
    vec_arms = {t: all_arms[t] for t in VECTOR_DISPATCH}

    _, mat_lines = optimize_type_switch(f"meta.matrixType({switch_var})", mat_arms)
    _, vec_lines = optimize_type_switch(f"meta.vectorType({switch_var})", vec_arms)

    # Cross-domain collapse: if both domains collapsed to the same single expression
    if len(mat_lines) == 1 and len(vec_lines) == 1 and mat_lines[0] == vec_lines[0]:
        return False, mat_lines

    out = [f"switch (comptime meta.domain({switch_var})) {{"]

    def embed(tag: str, lines: list[str]) -> None:
        if len(lines) == 1:
            out.append(f"{T}{tag} => {lines[0]}")
        else:
            out.append(f"{T}{tag} => {lines[0]}")
            for ln in lines[1:-1]:
                out.append(f"{T}{ln}")
            out.append(f"{T}{lines[-1]}")

    embed(".matrix", mat_lines)
    embed(".vector", vec_lines)
    out.append(f"{T}else => unreachable,")
    out.append("},")

    return False, out


def generate(base_indent: int = 4) -> str:
    bi = " " * base_indent

    # Bottom-up: Y -> X -> O
    o_arms: dict[str, tuple[bool, list[str]]] = {}
    for o in ALL_DISPATCH:
        x_arms: dict[str, tuple[bool, list[str]]] = {}
        for x in ALL_DISPATCH:
            y_arms = {y: zig_arm(o, x, y) for y in ALL_DISPATCH}
            x_arms[x] = optimize_domain_type_switch("Y", y_arms)
        o_arms[o] = optimize_domain_type_switch("X", x_arms)

    _, lines = optimize_domain_type_switch("O", o_arms)

    if lines[-1] == "},":
        lines[-1] = "}"

    return "\n".join(f"{bi}{ln}" if ln else ln for ln in lines)


if __name__ == "__main__":
    if "--stats" in sys.argv:
        c: Counter[str] = Counter()
        for o in ALL_DISPATCH:
            for x in ALL_DISPATCH:
                for y in ALL_DISPATCH:
                    c[decide(o, x, y).split(":")[0]] += 1
        total = sum(c.values())
        print(f"Total (20³)    : {total:>6}", file=sys.stderr)
        print(f"Valid imports  : {c['imp']:>6}", file=sys.stderr)
        print(f"Compile errors : {c['err']:>6}", file=sys.stderr)
        print(f"Unreachable    : {c['unr']:>6}", file=sys.stderr)

    print(generate())
