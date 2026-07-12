const std = @import("std");

const meta = @import("../meta.zig");

/// Determines whether a matrix type utilizes compressed pointer-and-index
/// storage (CSR or CSC).
pub fn isCompressed(comptime M: type) bool {
    return comptime switch (meta.matrixType(M)) {
        .general_sparse, .symmetric_sparse, .hermitian_sparse, .triangular_sparse => true,
        else => false,
    };
}

/// Evaluates whether a triangular matrix type implicitly represents a unit
/// diagonal.
pub fn hasImplicitDiag(comptime M: type) bool {
    return comptime meta.matrixType(M) == .triangular_sparse and meta.diagOf(M) == .unit;
}

/// Assesses whether a symmetric or Hermitian operand requires its implicit
/// off-diagonal mirror materialized to satisfy the structure of `O`.
pub fn needsMirror(comptime O: type, comptime M: type) bool {
    const mt = comptime meta.matrixType(M);
    if (comptime mt != .symmetric_sparse and mt != .hermitian_sparse)
        return false;

    return comptime meta.isSparseVector(O) or meta.matrixType(O) == .general_sparse;
}

/// Performs an O(log K) binary search within a compressed storage line to
/// verify the existence of an explicit coordinate.
pub fn hasExplicit(m: anytype, row: usize, col: usize) bool {
    const M = @TypeOf(m);

    const row_major = comptime meta.layoutOf(M) == .row_major;
    const line = if (row_major) row else col;
    const target = if (row_major) col else row;

    var lo: usize = m.ptr[line];
    var hi: usize = m.ptr[line + 1];
    while (lo < hi) {
        const mid = lo + (hi - lo) / 2;
        if (m.idx[mid] == target) return true;
        if (m.idx[mid] < target) lo = mid + 1 else hi = mid;
    }

    return false;
}

/// Evaluates whether a matrix instance contains a mathematical non-zero at
/// `(row, col)`, accounting for its specific storage layout and implicit
/// properties.
pub fn hasNZ(m: anytype, row: usize, col: usize) bool {
    const M = @TypeOf(m);

    return switch (comptime meta.matrixType(M)) {
        .diagonal_static, .diagonal_sparse => row == col,
        .permutation_static, .permutation_sparse => switch (comptime M.direction) {
            .forward => m.idx[row] == col,
            .backward => m.idx[col] == row,
        },
        .general_sparse => hasExplicit(m, row, col),
        .symmetric_sparse, .hermitian_sparse => hasExplicit(m, row, col) or hasExplicit(m, col, row),
        .triangular_sparse => ((comptime meta.diagOf(M) == .unit) and row == col) or hasExplicit(m, row, col),
        else => unreachable,
    };
}

pub fn format(formatter: anytype, comptime num_fmt: []const u8, rows: usize, cols: usize, writer: *std.Io.Writer) !void {
    const gap_size: usize = 3;
    const gap_str = " " ** gap_size;

    var cell_width: usize = 0;
    for (0..rows) |i| {
        for (0..cols) |j| {
            const val = formatter.matrix.get(i, j) catch unreachable;
            const val_len = std.fmt.count(num_fmt, .{val});
            if (val_len > cell_width)
                cell_width = val_len;
        }
    }

    if (cell_width == 0)
        cell_width = 1;

    const inner_width = (cols * cell_width) + ((cols - 1) * gap_size);

    // Top border
    try writer.writeAll("          ┌");
    for (0..inner_width + 2) |_| try writer.writeAll(" ");
    try writer.writeAll("┐\n");

    // Data
    for (0..rows) |i| {
        try writer.writeAll("          │ ");
        for (0..cols) |j| {
            if (j > 0) try writer.writeAll(gap_str);
            const val = formatter.matrix.get(i, j) catch unreachable;
            const val_len = std.fmt.count(num_fmt, .{val});

            for (0..(cell_width - val_len)) |_| try writer.writeAll(" ");
            try writer.print(num_fmt, .{val});
        }
        try writer.writeAll(" │\n");
    }

    // Bot Border
    try writer.writeAll("          └");
    for (0..inner_width + 2) |_| try writer.writeAll(" ");
    try writer.writeAll("┘");
}
