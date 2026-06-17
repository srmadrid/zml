const std = @import("std");

const zsl = @import("zsl");

const tzsl = @import("zsl.zig");

pub const general = struct {
    pub fn Static(N: type, layout: zsl.matrix.Layout) type {
        return struct {
            pub const instantiate = true;
            pub const general = true;
            pub const storage_layout = layout;
            pub const Numeric = N;
        };
    }
};

pub const symmetric = struct {
    pub fn Static(N: type, uplo: zsl.matrix.Uplo, layout: zsl.matrix.Layout) type {
        return struct {
            pub const instantiate = true;
            pub const symmetric = true;
            pub const storage_uplo = uplo;
            pub const storage_layout = layout;
            pub const Numeric = N;
        };
    }
};

pub const hermitian = struct {
    pub fn Static(N: type, uplo: zsl.matrix.Uplo, layout: zsl.matrix.Layout) type {
        return struct {
            pub const instantiate = true;
            pub const hermitian = true;
            pub const storage_uplo = uplo;
            pub const storage_layout = layout;
            pub const Numeric = N;
        };
    }
};

pub const triangular = struct {
    pub fn Static(N: type, uplo: zsl.matrix.Uplo, diag: zsl.matrix.Diag, layout: zsl.matrix.Layout) type {
        return struct {
            pub const instantiate = true;
            pub const triangular = true;
            pub const storage_uplo = uplo;
            pub const storage_diag = diag;
            pub const storage_layout = layout;
            pub const Numeric = N;
        };
    }
};

pub const diagonal = struct {
    pub fn Static(N: type) type {
        return struct {
            pub const instantiate = true;
            pub const diagonal = true;
            pub const Numeric = N;
        };
    }
};

pub const permutation = struct {
    pub fn Static(N: type, direction: zsl.matrix.permutation.Direction) type {
        return struct {
            pub const instantiate = true;
            pub const permutation = true;
            pub const storage_direction = direction;
            pub const Numeric = N;
        };
    }
};

pub fn printMatrix(desc: []const u8, a: anytype) void {
    const A = @TypeOf(a);

    std.debug.print("\nMatrix {s}:\n\n", .{desc});

    const rows = if (comptime zsl.meta.isStaticMatrix(A)) A.rows else a.rows;
    const cols = if (comptime zsl.meta.isStaticMatrix(A)) A.cols else a.cols;

    var i: usize = 0;
    while (i < rows) : (i += 1) {
        std.debug.print("\t", .{});

        var j: usize = 0;
        while (j < cols) : (j += 1) {
            if (comptime zsl.meta.isComplex(zsl.meta.Numeric(@TypeOf(a)))) {
                std.debug.print("{d} + {d}i    ", .{ (a.get(i, j) catch unreachable).re, (a.get(i, j) catch unreachable).im });
            } else {
                std.debug.print("{d}    ", .{a.get(i, j) catch unreachable});
            }
        }
        std.debug.print("\n", .{});
    }
    std.debug.print("\n", .{});
}

fn randomPermutation(rand: std.Random, data: []usize) void {
    // Initialize with identity permutation
    var i: usize = 0;
    while (i < data.len) : (i += 1) {
        data[i] = i;
    }

    // Shuffle using Fisher-Yates algorithm
    i = data.len - 1;
    while (i > 0) : (i -= 1) {
        const j = rand.intRangeAtMost(usize, 0, i);
        const temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }
}

pub fn randomMatrix(comptime M: type, allocator: std.mem.Allocator, rand: std.Random, comptime rows: usize, comptime cols: usize, nnz: usize) !M {
    switch (comptime zsl.meta.matrixType(M)) {
        .general_static => {
            var result: M = .init;

            inline for (0..rows) |i| {
                inline for (0..cols) |j| {
                    result.set(
                        i,
                        j,
                        tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                    ) catch unreachable;
                }
            }

            return result;
        },
        .general_dense => {
            var result: M = try .init(allocator, rows, cols);

            var i: usize = 0;
            while (i < rows) : (i += 1) {
                var j: usize = 0;
                while (j < cols) : (j += 1) {
                    result.set(
                        i,
                        j,
                        tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                    ) catch unreachable;
                }
            }

            return result;
        },
        .general_sparse => {
            if (nnz == rows * cols)
                return M.init(allocator, rows, cols, nnz);

            var builder: zsl.matrix.builder.Sparse(zsl.meta.Numeric(M)) = try .init(allocator, rows, cols, nnz);
            errdefer builder.deinit(allocator);

            var count: usize = 0;
            while (count < nnz) : (count += 1) {
                builder.appendAssumeCapacity(
                    rand.intRangeAtMost(usize, 0, rows - 1),
                    rand.intRangeAtMost(usize, 0, cols - 1),
                    tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                );
            }

            return builder.compile(allocator, M);
        },
        .symmetric_static, .hermitian_static => {
            var result: M = .init;

            inline for (0..rows) |i| {
                inline for (i..cols) |j| {
                    result.set(
                        i,
                        j,
                        zsl.numeric.cast(
                            zsl.meta.Numeric(M),
                            if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                                zsl.cf64{ .re = rand.float(f64), .im = if ((comptime zsl.meta.isHermitianMatrix(M)) and i == j) 0.0 else rand.float(f64) }
                            else
                                rand.float(f64),
                        ),
                    ) catch unreachable;
                }
            }

            return result;
        },
        .symmetric_dense, .hermitian_dense => {
            var result: M = try .init(allocator, rows);

            var i: usize = 0;
            while (i < rows) : (i += 1) {
                var j: usize = i;
                while (j < rows) : (j += 1) {
                    result.set(
                        i,
                        j,
                        zsl.numeric.cast(
                            zsl.meta.Numeric(M),
                            if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                                zsl.cf64{ .re = rand.float(f64), .im = if ((comptime zsl.meta.isHermitianMatrix(M)) and i == j) 0.0 else rand.float(f64) }
                            else
                                rand.float(f64),
                        ),
                    ) catch unreachable;
                }
            }

            return result;
        },
        .symmetric_sparse, .hermitian_sparse => {
            if (nnz == rows * cols)
                return M.init(allocator, rows, nnz);

            var builder: zsl.matrix.builder.Sparse(zsl.meta.Numeric(M)) = try .init(allocator, rows, rows, nnz);
            errdefer builder.deinit(allocator);

            var count: usize = 0;
            while (count < nnz) : (count += 1) {
                const r = rand.intRangeAtMost(usize, 0, rows - 1);
                const c = rand.intRangeAtMost(usize, 0, cols - 1);

                builder.appendAssumeCapacity(
                    r,
                    c,
                    zsl.numeric.cast(
                        zsl.meta.Numeric(M),
                        if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                            zsl.cf64{ .re = rand.float(f64), .im = if ((comptime zsl.meta.isHermitianMatrix(M)) and r == c) 0.0 else rand.float(f64) }
                        else
                            rand.float(f64),
                    ),
                );
            }

            return builder.compile(allocator, M);
        },
        .triangular_static => {
            var result: M = .init;

            if (comptime zsl.meta.uploOf(M) == .upper) {
                inline for (0..(comptime zsl.int.min(rows, cols))) |i| {
                    if (comptime zsl.meta.diagOf(M) == .non_unit) {
                        result.set(
                            i,
                            i,
                            tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                        ) catch unreachable;
                    }

                    inline for (i + 1..cols) |j| {
                        result.set(
                            i,
                            j,
                            tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                        ) catch unreachable;
                    }
                }
            } else {
                inline for (0..rows) |i| {
                    inline for (0..(comptime zsl.int.min(i, cols))) |j| {
                        result.set(
                            i,
                            j,
                            tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                        ) catch unreachable;
                    }

                    if ((comptime zsl.meta.diagOf(M) == .non_unit) and i < cols) {
                        result.set(
                            i,
                            i,
                            tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                        ) catch unreachable;
                    }
                }
            }

            return result;
        },
        .triangular_dense => {
            var result: M = try M.init(allocator, rows, cols);

            if (comptime zsl.meta.uploOf(M) == .upper) {
                var i: usize = 0;
                while (i < zsl.int.min(rows, cols)) : (i += 1) {
                    if (comptime zsl.meta.diagOf(M) == .non_unit) {
                        result.set(
                            i,
                            i,
                            tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                        ) catch unreachable;
                    }

                    var j: usize = i + 1;
                    while (j < cols) : (j += 1) {
                        result.set(
                            i,
                            j,
                            tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                        ) catch unreachable;
                    }
                }
            } else {
                var i: usize = 0;
                while (i < rows) : (i += 1) {
                    var j: usize = 0;
                    while (j < i and j < cols) : (j += 1) {
                        result.set(
                            i,
                            j,
                            tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                        ) catch unreachable;
                    }

                    if ((comptime zsl.meta.diagOf(M) == .non_unit) and i < cols) {
                        result.set(
                            i,
                            i,
                            tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                        ) catch unreachable;
                    }
                }
            }

            return result;
        },
        .triangular_sparse => {
            if (nnz == rows * cols)
                return M.init(allocator, rows, cols, nnz);

            var builder: zsl.matrix.builder.Sparse(zsl.meta.Numeric(M)) = try .init(allocator, rows, cols, nnz);
            errdefer builder.deinit(allocator);

            var count: usize = 0;
            while (count < nnz) : (count += 1) {
                const r = rand.intRangeAtMost(usize, 0, rows - 1);
                const c = rand.intRangeAtMost(usize, 0, cols - 1);

                builder.appendAssumeCapacity(
                    r,
                    c,
                    tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                );
            }

            return builder.compile(allocator, M);
        },
        .diagonal_static => {
            var result: M = .init;

            inline for (0..(comptime zsl.int.min(rows, cols))) |i| {
                result.set(
                    i,
                    i,
                    tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                ) catch unreachable;
            }

            return result;
        },
        .diagonal_sparse => {
            var result: M = try .init(allocator, rows, cols);
            errdefer result.deinit(allocator);

            var i: usize = 0;
            while (i < zsl.int.min(rows, cols)) : (i += 1) {
                result.set(
                    i,
                    i,
                    tzsl.randomNumber(zsl.meta.Numeric(M), rand),
                ) catch unreachable;
            }

            return result;
        },
        .builder_sparse => return .init(allocator, rows, cols, nnz),
        .permutation_static => {
            var result: M = .init;

            randomPermutation(rand, result.data[0..rows]);

            return result;
        },
        .permutation_sparse => {
            var result: M = try .init(allocator, rows);
            errdefer result.deinit(allocator);

            randomPermutation(rand, result.data[0..rows]);

            return result;
        },
        .numeric => unreachable,
    }
}

pub fn correctApply2(comptime O: type, allocator: std.mem.Allocator, rows: usize, cols: usize, A: anytype, B: anytype, op_: anytype) !zsl.matrix.general.Dense(O, .col_major) {
    const result: zsl.matrix.general.Dense(O, .col_major) = try .init(allocator, rows, cols);

    var j: usize = 0;
    while (j < result.cols) : (j += 1) {
        var i: usize = 0;
        while (i < result.rows) : (i += 1) {
            op_(
                &result.data[result._index(i, j)],
                if (comptime zsl.meta.isMatrix(@TypeOf(A))) A.get(i, j) catch unreachable else A,
                if (comptime zsl.meta.isMatrix(@TypeOf(B))) B.get(i, j) catch unreachable else B,
            );
        }
    }

    return result;
}

pub fn areEql(A: anytype, B: anytype) !void {
    var all_eql = true;

    var j: usize = 0;
    while (j < B.cols) : (j += 1) {
        var i: usize = 0;
        while (i < B.rows) : (i += 1) {
            all_eql = all_eql and zsl.numeric.eq(A.get(i, j) catch unreachable, B.get(i, j) catch unreachable);
        }
    }

    if (!all_eql)
        return error.NotEqual;
}

test {
    _ = @import("matrix/apply2.zig");
    _ = @import("matrix/apply2Alloc.zig");
    _ = @import("matrix/apply2Into.zig");
}
