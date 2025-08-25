const std = @import("std");
const zml = @import("zml");
const ci = @cImport({
    @cInclude("lapacke.h");
});

fn avg(values: []const f64) f64 {
    var sum: f64 = 0;
    for (values) |value| {
        sum += value;
    }
    return zml.float.div(sum, values.len);
}

fn avg_complex(values: []const zml.cf64) f64 {
    var sum: f64 = 0;
    for (values) |value| {
        sum += zml.cfloat.abs(value);
    }
    return zml.float.div(sum, values.len);
}

fn random_buffer(
    allocator: std.mem.Allocator,
    size: u32,
) ![]f64 {
    var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
    const rand = prng.random();

    const buffer = try allocator.alloc(f64, size);
    for (0..buffer.len) |i| {
        buffer[i] = rand.float(f64);
    }
    return buffer;
}

/// Generate a random m×n matrix A with specified 2-norm condition number `kappa`.
pub fn random_matrix(
    allocator: std.mem.Allocator,
    m: u32,
    n: u32,
    kappa: f64,
    order: zml.Order,
) ![]f64 {
    // - `allocator`: memory allocator
    // - `m`, `n`: dimensions
    // - `kappa`: desired condition number κ₂(A)
    // - `order`: either .RowMajor or .ColMajor
    //
    // Method:
    // 1. Construct singular values σᵢ geometrically spaced between 1 and κ^(1/min(m,n)).
    // 2. Generate random m×min(m,n) matrix X, QR factor → Q₁ (m×r).
    // 3. Generate random n×min(m,n) matrix Y, QR factor → Q₂ (n×r).
    // 4. Form A = Q₁ * diag(σ) * Q₂ᵀ.
    //     - Compute B = diag(σ) * Q₂ᵀ.
    //     - Then A = Q₁ × B.
    const r = zml.int.min(m, n);

    // 1) geometric singular values σ[i]
    var sig = try allocator.alloc(f64, r);
    defer allocator.free(sig);
    for (0..r) |i| {
        sig[i] = zml.float.pow(kappa, zml.float.div(zml.scast(f64, i), r - 1));
    }

    // 2) random X: m×r → QR → Q₁
    const X = try random_buffer(allocator, m * r);
    defer allocator.free(X);
    const tauX = try allocator.alloc(f64, r);
    defer allocator.free(tauX);
    _ = ci.LAPACKE_dgeqrf(
        order.toCInt(),
        zml.scast(c_int, m),
        zml.scast(c_int, r),
        X.ptr,
        zml.scast(c_int, if (order == .col_major) m else r),
        tauX.ptr,
    );
    _ = ci.LAPACKE_dorgqr(
        order.toCInt(),
        zml.scast(c_int, m),
        zml.scast(c_int, r),
        zml.scast(c_int, r),
        X.ptr,
        zml.scast(c_int, if (order == .col_major) m else r),
        tauX.ptr,
    );
    // X now holds Q₁

    // 3) random Y: n×r → QR → Q₂
    const Y = try random_buffer(allocator, n * r);
    defer allocator.free(Y);
    const tauY = try allocator.alloc(f64, r);
    defer allocator.free(tauY);
    _ = ci.LAPACKE_dgeqrf(
        order.toCInt(),
        zml.scast(c_int, n),
        zml.scast(c_int, r),
        Y.ptr,
        zml.scast(c_int, if (order == .col_major) n else r),
        tauY.ptr,
    );
    _ = ci.LAPACKE_dorgqr(
        order.toCInt(),
        zml.scast(c_int, n),
        zml.scast(c_int, r),
        zml.scast(c_int, r),
        Y.ptr,
        zml.scast(c_int, if (order == .col_major) n else r),
        tauY.ptr,
    );
    // Y now holds Q₂

    // 4a) form B = diag(sig) * Q2ᵀ (B is r×n in same layout)
    const B = try allocator.alloc(f64, r * n);
    defer allocator.free(B);
    for (0..r) |i| {
        const row_start = if (order == .col_major)
            B.ptr + i
        else
            B.ptr + i * n;

        const stride = if (order == .col_major) r else 1;

        for (0..n) |j| {
            // Q2[j,i] is at Y.ptr + j*ldY + i
            const q2ji = (Y.ptr + j * if (order == .col_major) n else r)[i];
            row_start[j * stride] = sig[i] * q2ji;
        }
    }

    // 4b) A = Q1 (m×r) × B (r×n) → A (m×n)
    const A = try allocator.alloc(f64, m * n);
    zml.linalg.blas.dgemm(
        order,
        .no_trans,
        .no_trans,
        zml.scast(i32, m),
        zml.scast(i32, n),
        zml.scast(i32, r),
        1,
        X.ptr,
        zml.scast(i32, if (order == .col_major) m else r),
        B.ptr,
        zml.scast(i32, if (order == .col_major) r else n),
        0,
        A.ptr,
        zml.scast(i32, if (order == .col_major) m else n),
    );

    return A;
}

fn max_difference(a: []const f64, b: []const f64) struct {
    index: u32,
    value: f64,
} {
    std.debug.assert(a.len == b.len);

    var max_diff: f64 = 0;
    var max_index: u32 = 0;

    var i: u32 = 0;
    while (i < a.len) : (i += 1) {
        const diff = zml.float.abs(a[i] - b[i]);
        if (diff > max_diff) {
            max_diff = diff;
            max_index = i;
        }
    }

    return .{ .index = max_index, .value = max_diff };
}

fn random_symmetric_matrix(
    allocator: std.mem.Allocator,
    size: u32,
    factor: f64,
) ![]f64 {
    var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
    const rand = prng.random();

    const matrix = try allocator.alloc(f64, size * size);

    for (0..size) |i| {
        for (i + 1..size) |j| {
            const value = rand.float(f64) * factor;
            matrix[i * size + j] = value;
            matrix[j * size + i] = value;
        }
    }

    return matrix;
}

/// Generates a random symmetric positive definite matrix with a specified condition number.
fn random_symmetric_positive_definite_matrix(
    allocator: std.mem.Allocator,
    size: u32,
    cond_target: f64,
) ![]f64 {
    // 1) compute λ_i geometrically between 1 and cond_target
    var lambdas = try allocator.alloc(f64, size);
    defer allocator.free(lambdas);
    const l_min = 1;
    const l_max = cond_target;
    for (0..size) |i| {
        lambdas[i] = l_min * zml.float.pow(l_max / l_min, zml.float.div(zml.scast(f64, i), size - 1));
    }

    // 2) random X, then QR -> get Q (n×n)
    const X = try random_buffer(allocator, size * size);
    defer allocator.free(X);
    const tau = try allocator.alloc(f64, size);
    defer allocator.free(tau);
    // QR factorization in-place: X = Q*R
    _ = ci.LAPACKE_dgeqrf(
        zml.Order.col_major.toCInt(),
        zml.scast(c_int, size),
        zml.scast(c_int, size),
        X.ptr,
        zml.scast(c_int, size),
        tau.ptr,
    );
    _ = ci.LAPACKE_dorgqr(
        zml.Order.col_major.toCInt(),
        zml.scast(c_int, size),
        zml.scast(c_int, size),
        zml.scast(c_int, size),
        X.ptr,
        zml.scast(c_int, size),
        tau.ptr,
    );
    // now X holds Q

    // 3) form B = D * Q^T, where D = diag(λ_1, λ_2, ..., λ_n)
    var B = try allocator.alloc(f64, size * size);
    defer allocator.free(B);
    // copy Q^T into B
    for (0..size) |i| {
        for (0..size) |j|
            B[i * size + j] = X[j * size + i];
    }
    // scale row i of B by λ_i
    for (0..size) |i| {
        zml.linalg.blas.dscal(zml.scast(i32, size), lambdas[i], B.ptr + i * size, 1);
    }

    // 4) Assemble A = Q * D * Q^T as A = Q * B
    const A = try allocator.alloc(f64, size * size);
    zml.linalg.blas.dgemm(
        .col_major,
        .no_trans,
        .no_trans,
        zml.scast(i32, size),
        zml.scast(i32, size),
        zml.scast(i32, size),
        1,
        X.ptr,
        zml.scast(i32, size),
        B.ptr,
        zml.scast(i32, size),
        0,
        A.ptr,
        zml.scast(i32, size),
    );

    return A;
}

fn frobernius_norm_difference(a: []const f64, b: []const f64) f64 {
    std.debug.assert(a.len == b.len);

    var norm: f64 = 0;
    for (0..a.len) |i| {
        const diff = a[i] - b[i];
        norm += diff * diff;
    }
    return zml.float.sqrt(norm);
}

fn is_symmetric(a: []const f64, size: u32) bool {
    for (0..size) |i| {
        for (i + 1..size) |j| {
            if (!std.math.approxEqRel(f64, a[i * size + j], a[j * size + i], 1e-9)) {
                return false;
            }
        }
    }
    return true;
}

fn random_complex_matrix(
    allocator: std.mem.Allocator,
    rows: u32,
    cols: u32,
    factor: f64,
) ![]zml.cf64 {
    var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
    const rand = prng.random();

    const matrix = try allocator.alloc(zml.cf64, rows * cols);
    for (0..matrix.len) |i| {
        matrix[i] = zml.cf64.init(rand.float(f64) * factor, rand.float(f64) * factor);
    }
    return matrix;
}

fn random_complex_hermitian_positive_definite_matrix(
    allocator: std.mem.Allocator,
    size: u32,
    factor: f64,
) ![]zml.cf64 {
    // allocate M
    const M = try random_complex_matrix(allocator, size, size, factor);
    // allocate A = M M^H
    const A = try allocator.alloc(zml.cf64, size * size);

    // A = M × M^H
    zml.linalg.blas.zgemm(
        .row_major,
        .no_trans,
        .conj_trans,
        zml.scast(i32, size),
        zml.scast(i32, size),
        zml.scast(i32, size),
        zml.cf64.init(1, 0),
        M.ptr,
        zml.scast(i32, size),
        M.ptr,
        zml.scast(i32, size),
        zml.cf64.init(0, 0),
        A.ptr,
        zml.scast(i32, size),
    );

    allocator.free(M);
    return A;
}

fn frobenius_norm_complex_difference(a: []const zml.cf64, b: []const zml.cf64) f64 {
    std.debug.assert(a.len == b.len);

    var norm: f64 = 0;
    for (0..a.len) |i| {
        const diff = zml.cfloat.sub(a[i], b[i]);
        norm += diff.re * diff.re + diff.im * diff.im;
    }
    return zml.float.sqrt(norm);
}

fn is_hermitian(a: []const zml.cf64, size: u32) bool {
    for (0..size) |i| {
        for (i + 1..size) |j| {
            if (!std.math.approxEqRel(f64, a[i * size + j].re, a[j * size + i].re, 1e-9) or
                !std.math.approxEqRel(f64, a[i * size + j].im, -a[j * size + i].im, 1e-9))
            {
                return false;
            }
        }
    }
    return true;
}

fn print_matrix(desc: []const u8, m: u32, n: u32, a: []f64, lda: u32, order: zml.Order) void {
    std.debug.print("\n{s}\n", .{desc});
    if (order == .row_major) {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                std.debug.print("{d}  ", .{a[i * lda + j]});
            }
            std.debug.print("\n", .{});
        }
    } else {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                std.debug.print("{d}  ", .{a[i + j * lda]});
            }
            std.debug.print("\n", .{});
        }
    }
}

fn print_complex_matrix(desc: []const u8, m: u32, n: u32, a: []zml.cf64, lda: u32, order: zml.Order) void {
    std.debug.print("\n{s}\n", .{desc});
    if (order == .row_major) {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                std.debug.print("{d:.4} + {d:.4}i  ", .{ a[i * lda + j].re, a[i * lda + j].im });
            }
            std.debug.print("\n", .{});
        }
    } else {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                std.debug.print("{d:.4} + {d:.4}i  ", .{ a[i + j * lda].re, a[i + j * lda].im });
            }
            std.debug.print("\n", .{});
        }
    }
}

pub fn main() !void {
    // const a: std.mem.Allocator = std.heap.page_allocator;
    var gpa = std.heap.DebugAllocator(.{}){};
    defer _ = gpa.deinit();
    const a = gpa.allocator();
    //_ = a;

    // const a: u64 = 1000;
    // const b: f64 = 1000;
    // const c = zml.add(a, b, .{}) catch unreachable;

    // std.debug.print("{d} + {d} = {d}\n", .{ a, b, c });
    // std.debug.print("@TypeOf(c) = {s}\n", .{@typeName(@TypeOf(c))});

    // const s = zml.types.mixStructs(.{ .a = a, .b = b }, .{ .c = c });
    // std.debug.print("Mixed struct: {}\n", .{s});
    // std.debug.print("Type of s: {s}\n", .{@typeName(@TypeOf(s))});
    // std.debug.print("Fields of s: \n", .{});
    // inline for (@typeInfo(@TypeOf(s)).@"struct".fields) |field| {
    //     std.debug.print("\t{s}: {s}\n", .{ field.name, @typeName(field.type) });
    // }

    // const t = zml.types.stripStruct(s, &.{ "a", "c" });
    // std.debug.print("Stripped struct: {}\n", .{t});
    // std.debug.print("Type of t: {s}\n", .{@typeName(@TypeOf(t))});
    // std.debug.print("Fields of t: \n", .{});
    // inline for (@typeInfo(@TypeOf(t)).@"struct".fields) |field| {
    //     std.debug.print("\t{s}: {s}\n", .{ field.name, @typeName(field.type) });
    // }

    // try symbolicTesting(a);

    // try generalTesting(a);

    // try addTesting(a);

    // try iterTesting(a);

    // try iterPerfTesting(a);

    // try multiIterTesting(a);

    try perfTesting(a);

    // try typeTesting(a);

    // try transposeTesting(a);

    // try blasPerfTesting(a);

    // try lapackTesting(a);

    // try lapackPerfTesting(a);

    // try matrixTesting(a);

    // coreTesting();
}

fn ask_user(default: u32) !u32 {
    const stdin = std.io.getStdIn().reader();
    const stdout = std.io.getStdOut().writer();

    var buf: [10]u8 = undefined;

    try stdout.print("Enter the number of iterations: ", .{});

    if (try (stdin.readUntilDelimiterOrEof(buf[0..], '\n'))) |user_input| {
        return std.fmt.parseInt(u32, user_input, 10);
    } else {
        return default;
    }
}

fn fillMatrix(a: anytype, factor: u32) void {
    const A: type = zml.types.Numeric(@TypeOf(a));
    switch (comptime zml.types.matrixType(@TypeOf(a))) {
        .general, .triangular => {
            var i: u32 = 0;
            while (i < a.rows * a.cols) : (i += 1) {
                if (comptime zml.types.isComplex(A)) {
                    a.data[i] = A.init(zml.scast(zml.types.Scalar(A), i + factor), zml.scast(zml.types.Scalar(A), i + factor));
                } else {
                    a.data[i] = zml.scast(A, i + factor);
                }
            }
        },
        .symmetric, .hermitian => {
            var i: u32 = 0;
            while (i < a.size * a.size) : (i += 1) {
                if (i % (a.size + 1) == 0) {
                    if (comptime zml.types.isComplex(A)) {
                        a.data[i] = A.init(zml.scast(zml.types.Scalar(A), i + factor), if (comptime zml.types.isHermitianMatrix(@TypeOf(a))) 0 else zml.scast(zml.types.Scalar(A), i + factor));
                    } else {
                        a.data[i] = zml.scast(A, i + factor);
                    }
                } else {
                    if (comptime zml.types.isComplex(A)) {
                        a.data[i] = A.init(zml.scast(zml.types.Scalar(A), i + factor), zml.scast(zml.types.Scalar(A), i + factor));
                    } else {
                        a.data[i] = zml.scast(A, i + factor);
                    }
                }
            }
        },
        .diagonal => {
            var i: u32 = 0;
            while (i < a.size) : (i += 1) {
                if (comptime zml.types.isComplex(A)) {
                    a.data[i] = A.init(zml.scast(zml.types.Scalar(A), i + factor), zml.scast(zml.types.Scalar(A), i + factor));
                } else {
                    a.data[i] = zml.scast(A, i + factor);
                }
            }
        },
        .banded => {
            var i: u32 = 0;
            while (i < (a.lower + a.upper + 1) * if (zml.types.orderOf(@TypeOf(a)) == .col_major) a.cols else a.rows) : (i += 1) {
                if (comptime zml.types.isComplex(A)) {
                    a.data[i] = A.init(zml.scast(zml.types.Scalar(A), i + factor), zml.scast(zml.types.Scalar(A), i + factor));
                } else {
                    a.data[i] = zml.scast(A, i + factor);
                }
            }
        },
        .tridiagonal => {
            var i: u32 = 0;
            while (i < (3 * a.size - 2)) : (i += 1) {
                if (comptime zml.types.isComplex(A)) {
                    a.data[i] = A.init(zml.scast(zml.types.Scalar(A), i + factor), zml.scast(zml.types.Scalar(A), i + factor));
                } else {
                    a.data[i] = zml.scast(A, i + factor);
                }
            }
        },
        else => unreachable,
    }
}

fn printMatrix(name: []const u8, a: anytype) void {
    std.debug.print("Matrix {s}:\n", .{name});
    if (comptime zml.types.isSymmetricMatrix(@TypeOf(a)) or zml.types.isHermitianMatrix(@TypeOf(a)) or zml.types.isTridiagonalMatrix(@TypeOf(a))) {
        var i: u32 = 0;
        while (i < a.size) : (i += 1) {
            var j: u32 = 0;
            while (j < a.size) : (j += 1) {
                if (comptime zml.types.isComplex(zml.types.Numeric(@TypeOf(a)))) {
                    std.debug.print("{d:3} + {d:3}i  ", .{ (a.get(i, j) catch unreachable).re, (a.get(i, j) catch unreachable).im });
                } else {
                    std.debug.print("{d:3}  ", .{a.get(i, j) catch unreachable});
                }
            }
            std.debug.print("\n", .{});
        }
        std.debug.print("\n", .{});
    } else {
        var i: u32 = 0;
        while (i < a.rows) : (i += 1) {
            var j: u32 = 0;
            while (j < a.cols) : (j += 1) {
                if (comptime zml.types.isComplex(zml.types.Numeric(@TypeOf(a)))) {
                    std.debug.print("{d:3} + {d:3}i  ", .{ (a.get(i, j) catch unreachable).re, (a.get(i, j) catch unreachable).im });
                } else {
                    std.debug.print("{d:3}  ", .{a.get(i, j) catch unreachable});
                }
            }
            std.debug.print("\n", .{});
        }
        std.debug.print("\n", .{});
    }
}

fn perfTesting(a: std.mem.Allocator) !void {
    var A: zml.matrix.Symmetric(zml.cf64, .lower, .row_major) = try .init(a, 5);
    defer A.deinit(a);

    fillMatrix(A, 1);
    printMatrix("A", A);

    var B: zml.matrix.Hermitian(zml.cf64, .lower, .row_major) = try .init(a, 5);
    defer B.deinit(a);

    fillMatrix(B, 2);
    printMatrix("B", B);

    const start_time = std.time.nanoTimestamp();
    var C: zml.matrix.General(zml.cf64, .row_major) = try zml.matrix.apply2(a, A, B, zml.add, .{});
    const end_time: i128 = std.time.nanoTimestamp();
    defer C.deinit(a);

    std.debug.print("Took: {d} seconds\n\n", .{zml.float.div(end_time - start_time, 1e9)});

    printMatrix("C = A + B", C);
}

fn matrixTesting(a: std.mem.Allocator) !void {
    var A: zml.array.Dense(f64) = try .init(a, &.{ 10, 5, 8, 9 }, .{ .order = .col_major });
    defer A.deinit(a);

    var i: u32 = 0;
    while (i < A.size) : (i += 1) {
        A.data[i] = zml.scast(f64, i + 1);
        //A.data[i] = zml.cf64.init(zml.scast(f64, i + 1), zml.scast(f64, i + 1));
    }

    std.debug.print("Matrix A (as Dense):\n", .{});
    i = 0;
    while (i < A.shape[0]) : (i += 1) {
        var j: u32 = 0;
        while (j < A.shape[1]) : (j += 1) {
            std.debug.print("{d:2}  ", .{try A.get(&.{ i, j, 1 })});
            //std.debug.print("{d:3} + {d:3}i  ", .{ (try A.get(i, j)).re, (try A.get(i, j)).im });
        }

        std.debug.print("\n", .{});
    }

    var AG: zml.matrix.General(f64) = try A.asGeneralMatrix(null, &.{ 1, 0 });
    defer AG.deinit(a);

    std.debug.print("Matrix A:\n", .{});
    i = 0;
    while (i < AG.rows) : (i += 1) {
        var j: u32 = 0;
        while (j < AG.cols) : (j += 1) {
            std.debug.print("{d:2}  ", .{try AG.get(i, j)});
            //std.debug.print("{d:3} + {d:3}i  ", .{ (try AD.get(&.{ i, j })).re, (try AD.get(&.{ i, j })).im });
        }

        std.debug.print("\n", .{});
    }
}

fn lapackTesting(a: std.mem.Allocator) !void {
    // const m = 3;
    // const n = 5;
    // const lda_row = n;
    // const lda_col = m;

    // // Row-major storage
    // const a_row: []f64 = try a.alloc(f64, m * n);
    // defer a.free(a_row);

    // @memcpy(a_row, &[_]f64{ // Initialize with some values
    //     2, 6, 8, 5, 4, 3, 1,
    //     1, 4, 7, 2, 5, 6, 2,
    //     3, 5, 9, 6, 1, 7, 3,
    //     4, 8, 6, 3, 2, 8, 4,
    //     5, 2, 1, 4, 3, 9, 5,
    // });

    // // Column-major storage (transpose indexing layout)
    // const a_col: []f64 = try a.alloc(f64, m * n);
    // defer a.free(a_col);

    // var r: u32 = 0;
    // while (r < m) : (r += 1) {
    //     var c: u32 = 0;
    //     while (c < n) : (c += 1) {
    //         a_col[r + c * lda_col] = a_row[r * lda_row + c];
    //     }
    // }

    // // Pivot array (LAPACK style: 1-based indexing)
    // const ipiv: []const i32 = &.{3}; // Swap row1<->row3, row2<->row4
    // const k1 = 1;
    // const k2 = 2;
    // const incx = 1;

    // print_matrix("Original Row-major", m, n, a_row, lda_row, .row_major);
    // print_matrix("Original Col-major", m, n, a_col, lda_col, .col_major);

    // // Apply laswp in row-major
    // try zml.linalg.lapack.laswp(.row_major, n, a_row.ptr, lda_row, k1, k2, ipiv.ptr, incx);

    // // Apply laswp in col-major
    // try zml.linalg.lapack.laswp(.col_major, n, a_col.ptr, lda_col, k1, k2, ipiv.ptr, incx);

    // print_matrix("After LASWP (Row-major)", m, n, a_row, lda_row, .row_major);
    // print_matrix("After LASWP (Col-major)", m, n, a_col, lda_col, .col_major);

    const m = 5;
    const n = 7;

    const lda_col = m;
    const lda_row = n;
    var ipiv_col_array: [m]i32 = undefined;
    const ipiv_col: []i32 = &ipiv_col_array;
    var ipiv_row_array: [m]i32 = undefined;
    const ipiv_row: []i32 = &ipiv_row_array;

    // Column-major input
    const a_col: []f64 = try a.alloc(f64, m * n);
    defer a.free(a_col);

    @memcpy(a_col, &[_]f64{
        2, 1, 3, 4, 5,
        6, 4, 5, 8, 2,
        8, 7, 9, 6, 1,
        5, 2, 6, 3, 4,
        4, 5, 1, 2, 3,
        3, 6, 7, 8, 9,
        1, 2, 3, 4, 5,
    });

    // Copy to row-major layout
    const a_row: []f64 = try a.alloc(f64, m * n);
    defer a.free(a_row);

    var r: u32 = 0;
    while (r < m) : (r += 1) {
        var c: u32 = 0;
        while (c < n) : (c += 1) {
            a_row[r * lda_row + c] = a_col[r + c * lda_col];
        }
    }

    // Print original
    print_matrix("Original A (col-major)", m, n, a_col, lda_col, .col_major);
    print_matrix("Original A (row-major)", m, n, a_row, lda_row, .row_major);

    // Column-major dgetrf2
    const info_col: i32 = try zml.linalg.lapack.getrf2(.col_major, m, n, a_col.ptr, lda_col, ipiv_col.ptr, .{});
    // Row-major wrapper
    const info_row: i32 = try zml.linalg.lapack.getrf2(.row_major, m, n, a_row.ptr, lda_row, ipiv_row.ptr, .{});

    // Print ipiv (1-based)
    std.debug.print("\nipiv col-major: ", .{});
    for (0..m) |i|
        std.debug.print("{d}, ", .{ipiv_col[i]});
    std.debug.print("\nipiv row-major: ", .{});
    for (0..m) |i|
        std.debug.print("{d}, ", .{ipiv_row[i]});
    std.debug.print("\n", .{});

    // Print LU factors
    print_matrix("LU in A (col-major)", m, n, a_col, lda_col, .col_major);
    print_matrix("LU in A (row-major)", m, n, a_row, lda_row, .row_major);

    std.debug.print("\ninfo_col = {d}, info_row = {d}\n", .{ info_col, info_row });
}

fn lapackPerfTesting(a: std.mem.Allocator) !void {
    const iter = 1;
    // const m = 300;
    const n = 3000;
    const nrhs = 3000;

    // A
    const a_original: []f64 = try random_symmetric_positive_definite_matrix(a, n, 2);
    //const a_original: []zml.cf64 = try random_complex_hermitian_positive_definite_matrix(a, n, 0.1);
    defer a.free(a_original);
    const a_mean: f64 = avg(a_original);
    //const a_mean: f64 = avg_complex(a_zml);
    const w: []f64 = try a.alloc(f64, n);
    defer a.free(w);
    _ = ci.LAPACKE_dsyev(
        zml.Order.row_major.toCInt(),
        'N',
        'L',
        zml.scast(c_int, n),
        a_original.ptr,
        zml.scast(c_int, n),
        w.ptr,
    );
    const a_cond: f64 = w[n - 1] / w[0];

    const a_zml: []f64 = try a.alloc(f64, n * n);
    //const a_zml: []zml.cf64 = try a.alloc(zml.cf64, n * n);
    defer a.free(a_zml);
    @memcpy(a_zml, a_original);

    const a_lapacke: []f64 = try a.alloc(f64, n * n);
    //const a_lapacke: []zml.cf64 = try a.alloc(zml.cf64, n * n);
    defer a.free(a_lapacke);
    @memcpy(a_lapacke, a_original);

    //const lda_col = n;
    const lda_row = n;

    // B
    const b_original: []f64 = try random_matrix(a, n, nrhs, 2, .row_major);
    //const b_original: []zml.cf64 = try random_complex_matrix(a, n, nrhs, 10);
    defer a.free(b_original);
    const b_mean: f64 = avg(b_original);
    const S: []f64 = try a.alloc(f64, zml.int.min(n, nrhs));
    defer a.free(S);
    const work: []f64 = try a.alloc(f64, zml.int.min(n, nrhs));
    defer a.free(work);
    _ = ci.LAPACKE_dgesvd(
        zml.Order.row_major.toCInt(),
        'N',
        'N',
        zml.scast(c_int, n),
        zml.scast(c_int, nrhs),
        b_original.ptr,
        zml.scast(c_int, nrhs),
        S.ptr,
        null,
        zml.scast(c_int, n),
        null,
        zml.scast(c_int, nrhs),
        work.ptr,
    );
    const b_cond: f64 = S[0] / S[zml.int.min(n, nrhs) - 1];

    const b_zml: []f64 = try a.alloc(f64, n * nrhs);
    defer a.free(b_zml);
    @memcpy(b_zml, b_original);

    const b_lapacke: []f64 = try a.alloc(f64, n * nrhs);
    defer a.free(b_lapacke);
    @memcpy(b_lapacke, b_original);

    //const ldb_col = n;
    const ldb_row = nrhs;

    // Pivot array (LAPACK style: 1-based indexing)
    // const ipiv_zml: []i32 = try a.alloc(i32, n);
    // defer a.free(ipiv_zml);
    // const ipiv_lapacke: []i32 = try a.alloc(i32, n);
    // defer a.free(ipiv_lapacke);

    var info_zml: i32 = 0;
    var start_time: i128 = std.time.nanoTimestamp();
    for (0..iter) |_| {
        info_zml = try zml.linalg.lapack.posv(
            .row_major,
            .lower,
            n,
            nrhs,
            a_zml.ptr,
            lda_row,
            b_zml.ptr,
            ldb_row,
            .{},
        );
        // info_zml = try zml.linalg.lapack.potrf(
        //     .col_major,
        //     .lower,
        //     n,
        //     a_zml.ptr,
        //     lda_col,
        //     .{},
        // );

        // try zml.linalg.lapack.potrs(
        //     .col_major,
        //     .lower,
        //     n,
        //     nrhs,
        //     a_zml.ptr,
        //     lda_col,
        //     b_zml.ptr,
        //     ldb_col,
        //     .{},
        // );
    }
    var end_time: i128 = std.time.nanoTimestamp();

    std.debug.print("Zml potrf took: {d} seconds\n", .{zml.float.div(end_time - start_time, 1e9 * iter)});

    var info_lapacke: i32 = 0;
    start_time = std.time.nanoTimestamp();
    for (0..iter) |_| {
        info_lapacke = zml.scast(i32, ci.LAPACKE_dposv(
            zml.Order.row_major.toCInt(),
            'L',
            zml.scast(c_int, n),
            zml.scast(c_int, nrhs),
            a_lapacke.ptr,
            zml.scast(c_int, lda_row),
            b_lapacke.ptr,
            zml.scast(c_int, ldb_row),
        ));
        // info_lapacke = zml.scast(i32, ci.LAPACKE_dpotrf(
        //     @intFromEnum(zml.Order.col_major),
        //     'L',
        //     zml.scast(c_int, n),
        //     a_lapacke.ptr,
        //     zml.scast(c_int, lda_col),
        // ));

        // _ = ci.LAPACKE_dpotrs(
        //     @intFromEnum(zml.Order.col_major),
        //     'L',
        //     zml.scast(c_int, n),
        //     zml.scast(c_int, nrhs),
        //     a_lapacke.ptr,
        //     zml.scast(c_int, lda_col),
        //     b_lapacke.ptr,
        //     zml.scast(c_int, ldb_col),
        // );
    }
    end_time = std.time.nanoTimestamp();

    std.debug.print("Lapacke potrf took: {d} seconds\n", .{zml.float.div(end_time - start_time, 1e9 * iter)});

    std.debug.print("info_zml = {d}, info_lapacke = {d}\n", .{ info_zml, info_lapacke });

    // Frobenius norm of the difference
    const a_norm: f64 = frobernius_norm_difference(a_zml, a_lapacke);
    const b_norm: f64 = frobernius_norm_difference(b_zml, b_lapacke);
    //const norm: f64 = frobenius_norm_complex_difference(a_zml, a_lapacke);

    std.debug.print("Frobenius norm of the difference of a: {d}\n", .{a_norm});
    std.debug.print("Frobenius norm of the difference of b: {d}\n", .{b_norm});
    std.debug.print("Mean of the original matrix A: {d}\n", .{a_mean});
    std.debug.print("Mean of the original matrix B: {d}\n", .{b_mean});

    //print_matrix("A (zml)", n, n, a_zml, lda_col, .col_major);
    //print_matrix("A (lapacke)", n, n, a_lapacke, lda_col, .col_major);
    //print_matrix("B (zml)", n, nrhs, b_zml, ldb_row, .row_major);
    //print_matrix("B (lapacke)", n, nrhs, b_lapacke, ldb_row, .row_major);

    const b_max_diff = max_difference(b_zml, b_lapacke);
    std.debug.print("Max difference in B:\n B (zml): {d}\n B (lapacke): {d}\n", .{
        b_zml[b_max_diff.index],
        b_lapacke[b_max_diff.index],
    });
    std.debug.print("Max difference value: {d}\n", .{b_max_diff.value});

    std.debug.print("Condition number of A: {}\n", .{a_cond});
    std.debug.print("Condition number of B: {}\n", .{b_cond});
}

fn blasPerfTesting(a: std.mem.Allocator) !void {
    const alpha: f64 = 2;
    const beta: f64 = 3;

    const m = 1000;
    const n = 1500;
    const k = 2000;

    const A: []f64 = try a.alloc(f64, m * k);
    defer a.free(A);
    const B: []f64 = try a.alloc(f64, k * n);
    defer a.free(B);
    const C: []f64 = try a.alloc(f64, m * n);
    defer a.free(C);

    for (0..A.len) |i| {
        A[i] = @floatFromInt(i + 1);
    }

    for (0..B.len) |i| {
        B[i] = @floatFromInt(i + 1);
    }

    for (0..C.len) |i| {
        C[i] = 0;
    }

    const start_time: i128 = std.time.nanoTimestamp();
    for (0..10) |_| {
        zml.linalg.blas.gemm(.col_major, .no_trans, .no_trans, m, n, k, alpha, A.ptr, m, B.ptr, k, beta, C.ptr, m, .{}) catch unreachable;
    }
    const end_time: i128 = std.time.nanoTimestamp();

    std.debug.print("zml.linalg.blas.gemm took: {d} seconds\n", .{zml.float.div(end_time - start_time, 1e9)});
}

fn transposeTesting(a: std.mem.Allocator) !void {
    if (false) {
        var A: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 6, 6 }, .{ .order = .RowMajor });
        defer A.deinit();
        const At = try A.transpose(null);

        std.debug.print("A dimentions = {}\n", .{A.shape.len});

        std.debug.print("A.shape = [  ", .{});
        for (A.shape[0..A.ndim]) |dim| {
            std.debug.print("{}  ", .{dim});
        }
        std.debug.print("]\n", .{});

        std.debug.print("A.strides = [  ", .{});
        for (A.strides[0..A.ndim]) |stride| {
            std.debug.print("{}  ", .{stride});
        }
        std.debug.print("]\n", .{});

        std.debug.print("A.size = {}\n", .{A.size});

        std.debug.print("At dimentions = {}\n", .{At.shape.len});

        std.debug.print("At.shape = [  ", .{});
        for (At.shape[0..At.ndim]) |dim| {
            std.debug.print("{}  ", .{dim});
        }
        std.debug.print("]\n", .{});

        std.debug.print("At.strides = [  ", .{});
        for (At.strides[0..At.ndim]) |stride| {
            std.debug.print("{}  ", .{stride});
        }
        std.debug.print("]\n", .{});

        std.debug.print("At.size = {}\n", .{At.size});

        for (0..A.size) |i| {
            A.data[i] = @floatFromInt(i + 1);
        }
        std.debug.print("A =\n", .{});
        for (0..A.shape[0]) |i| {
            std.debug.print("\t", .{});
            for (0..A.shape[1]) |j| {
                std.debug.print("{!d:.2}  ", .{A.get(&[_]u32{ i, j })});
            }
            std.debug.print("\n", .{});
        }

        std.debug.print("At =\n", .{});
        for (0..At.shape[0]) |i| {
            std.debug.print("\t", .{});
            for (0..At.shape[1]) |j| {
                std.debug.print("{!d:.2}  ", .{At.get(&[_]u32{ i, j })});
            }
            std.debug.print("\n", .{});
        }

        var B: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 6, 6 }, .{ .order = .ColumnMajor });
        defer B.deinit();

        try B.add(A, try A.transpose(null));
        // or try B.add(A, At);

        std.debug.print("B =\n", .{});
        for (0..B.shape[0]) |i| {
            std.debug.print("\t", .{});
            for (0..B.shape[1]) |j| {
                std.debug.print("{!d:.2}  ", .{B.get(&[_]u32{ i, j })});
            }
            std.debug.print("\n", .{});
        }

        const C = At.flatten();

        std.debug.print("C.shape = [  ", .{});
        for (C.shape[0..C.ndim]) |dim| {
            std.debug.print("{}  ", .{dim});
        }
        std.debug.print("]\n", .{});

        std.debug.print("C.size = {}\n", .{C.size});

        std.debug.print("C =\n", .{});
        for (0..C.size) |i| {
            std.debug.print("{!d:.2}  ", .{C.data[i]});
        }
        std.debug.print("\n", .{});

        var D: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 6, 6 }, .{});
        defer D.deinit();
        D.setAll(1);
        var E: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 6, 6 }, .{});
        defer E.deinit();
        E.setAll(0);
    }

    var x = try zml.Array(f64).init(a, &.{ 5, 9, 4, 8 }, .{});

    std.debug.print("x.shape = [\n", .{});
    for (0..x.ndim) |i| {
        std.debug.print("x.shape[{}] = {}\n", .{ i, x.shape[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x.strides = [\n", .{});
    for (0..x.ndim) |i| {
        std.debug.print("x.strides[{}] = {}\n", .{ i, x.metadata.dense.strides[i] });
    }
    std.debug.print("]\n\n", .{});

    for (0..x.shape[0]) |i| {
        for (0..x.shape[1]) |j| {
            for (0..x.shape[2]) |k| {
                for (0..x.shape[3]) |l| {
                    try x.set(&[_]u32{ i, j, k, l }, @floatFromInt(i * 1000 + j * 100 + k * 10 + l));
                }
            }
        }
    }

    const x1 = try x.slice(&.{ try .init(1, 5, 1), try .init(2, 5, 1), .single(2) });

    std.debug.print("x1 = x[1:5, 2:5, 2]\n", .{});
    std.debug.print("x1.shape = [\n", .{});
    for (0..x1.ndim) |i| {
        std.debug.print("x1.shape[{}] = {}\n", .{ i, x1.shape[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x1.strides = [\n", .{});
    for (0..x1.ndim) |i| {
        std.debug.print("x1.strides[{}] = {}\n", .{ i, x1.metadata.strided.strides[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x1.offset = {}\n\n", .{x1.metadata.strided.offset});

    const x2 = try x.slice(&.{ .all, try .init(2, 9, 1) });

    std.debug.print("x2 = x[:, 2:10]\n", .{});
    std.debug.print("x2.shape = [\n", .{});
    for (0..x2.ndim) |i| {
        std.debug.print("x2.shape[{}] = {}\n", .{ i, x2.shape[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x2.strides = [\n", .{});
    for (0..x2.ndim) |i| {
        std.debug.print("x2.strides[{}] = {}\n", .{ i, x2.metadata.strided.strides[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x2.offset = {}\n\n", .{x2.metadata.strided.offset});

    const x3 = try x.slice(&.{ try .init(1, 5, 1), try .init(8, 4, -3), .all_reverse, .all });

    std.debug.print("x3 = x[1:5, 8:4:-3, :, :]\n", .{});
    std.debug.print("x3.shape = [\n", .{});
    for (0..x3.ndim) |i| {
        std.debug.print("x3.shape[{}] = {}\n", .{ i, x3.shape[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x3.strides = [\n", .{});
    for (0..x3.ndim) |i| {
        std.debug.print("x3.strides[{}] = {}\n", .{ i, x3.metadata.strided.strides[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x3.offset = {}\n\n", .{x3.metadata.strided.offset});

    std.debug.print("x3 =\n", .{});
    for (0..x3.shape[2]) |k| {
        std.debug.print("Slice {}:\n", .{k});
        for (0..x3.shape[0]) |i| {
            std.debug.print("\t", .{});
            for (0..x3.shape[1]) |j| {
                std.debug.print("{d}  ", .{(try x3.get(&.{ i, j, k, 0 })).*});
            }
            std.debug.print("\n", .{});
        }
        std.debug.print("\n", .{});
    }
}

fn typeTesting(a: std.mem.Allocator) !void {
    const Complex = std.math.Complex;

    //const z = Complex(f64).init(10, 10);
    //const w = Complex(f64).init(20, 20);
    //var u: Complex(f64) = undefined;
    //const _add = @import("core/core.zig").supported._add;
    //_add(&u, z, w);
    //std.debug.print("u: {}\n", .{u});

    var B: zml.Array(Complex(f64)) = try zml.Array(Complex(f64)).init(a, &.{ 1, 1 }, .{ .order = .ColumnMajor });
    defer B.deinit();
    for (0..B.size) |i| {
        B.data[i].re = @floatFromInt(i + 1);
        B.data[i].im = @floatFromInt((i + 1) * 2);
    }
    std.debug.print("B =\n", .{});
    for (0..B.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..B.shape[1]) |j| {
            const z = try B.get(&.{ i, j });
            std.debug.print("{d:.2} + {d:.2}i  ", .{ z.re, z.im });
        }
        std.debug.print("\n", .{});
    }

    var C: zml.Array(Complex(f64)) = try zml.Array(Complex(f64)).init(a, &.{ 5, 8 }, .{ .order = .RowMajor });
    defer C.deinit();
    for (0..C.size) |i| {
        C.data[i].re = @floatFromInt(i + 1);
        C.data[i].im = @floatFromInt((i + 1) * 2);
    }
    std.debug.print("\nC =\n", .{});
    for (0..C.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..C.shape[1]) |j| {
            const z = try C.get(&.{ i, j });
            std.debug.print("{d:.2} + {d:.2}i  ", .{ z.re, z.im });
        }
        std.debug.print("\n", .{});
    }

    var D: zml.Array(Complex(f64)) = try zml.Array(Complex(f64)).init(a, &.{ 5, 8 }, .{ .order = .ColumnMajor });
    defer D.deinit();
    try D.add(B, C);
    // try zml.Array(u64).add(&D, B, C);
    std.debug.print("\nD =\n", .{});
    for (0..D.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..D.shape[1]) |j| {
            const z = try D.get(&.{ i, j });
            std.debug.print("{d:.2} + {d:.2}i  ", .{ z.re, z.im });
        }
        std.debug.print("\n", .{});
    }
}

fn symbolicTesting(a: std.mem.Allocator) !void {
    // Does not work, just for design.
    //const x = try zml.Variable.init(a, "x", &zml.Set.RealNumbers);
    //const S = try zml.Set.init.builder(a, "S", x, &[_]zml.Expression{zml.Expression.init.fromString(a, "sin(x) = 0", &[_]*zml.Symbol{&x})});
    //_ = S;
    const expr = "S = \\{x\\in\\mathbb{R}\\mid x > 0, \\arcsin(x) = 2\\pi k, \\forall k\\in\\mathbb{N}\\}";
    const arr = try zml.Expression.tokenize(a, expr);
    std.debug.print("{s}\n\nTokenized:\n", .{expr});
    for (arr.items) |token| {
        std.debug.print("<string = {s}, type = {}>\n", .{ token.string, token.type });
    }
}

fn generalTesting(a: std.mem.Allocator) !void {
    std.debug.print("Size of flags: {}\n", .{@sizeOf(zml.array.Flags)});

    var A: zml.Array(f64) = try .init(a, &.{ 20, 15, 8, 18 }, .{});
    defer A.deinit(a);

    std.debug.print("A dimentions = {}\n", .{A.ndim});

    std.debug.print("A.shape = [  ", .{});
    for (A.shape[0..A.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});

    std.debug.print("A.strides = [  ", .{});
    for (A.metadata.dense.strides[0..A.ndim]) |stride| {
        std.debug.print("{}  ", .{stride});
    }
    std.debug.print("]\n", .{});

    std.debug.print("A.size = {}\n", .{A.size});

    var B: zml.Array(f64) = try .logspace(a, 2, 3, 4, 10, false, .{});
    defer B.deinit(a);

    std.debug.print("B = [", .{});
    for (0..B.shape[0]) |i| {
        std.debug.print("{d}, ", .{B.data[i]});
    }
    std.debug.print("]\n", .{});

    var C: zml.Array(f64) = try zml.abs(B, .{ .array_allocator = a });
    defer C.deinit(a);

    std.debug.print("C = [", .{});
    for (0..C.shape[0]) |i| {
        std.debug.print("{d}, ", .{C.data[i]});
    }
    std.debug.print("]\n", .{});

    var D: zml.Array(zml.cf64) = try zml.Array(zml.cf64).init(a, &.{5}, .{});
    defer D.deinit(a);
    for (0..D.size) |i| {
        D.data[i] = zml.cf64.init(@floatFromInt(i + 1), @floatFromInt((i + 1) * 2));
    }

    std.debug.print("D = [", .{});
    for (0..D.shape[0]) |i| {
        std.debug.print("{d} + {d}i, ", .{ D.data[i].re, D.data[i].im });
    }
    std.debug.print("]\n", .{});

    const D_broadcasted: zml.Array(zml.cf64) = try D.broadcast(&.{ 3, 5 });

    std.debug.print("D_broadcasted.shape = [  ", .{});
    for (D_broadcasted.shape[0..D_broadcasted.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});
    std.debug.print("D_broadcasted = [\n", .{});
    for (0..D_broadcasted.shape[0]) |i| {
        for (0..D_broadcasted.shape[1]) |j| {
            const z = try D_broadcasted.get(&.{ i, j });
            std.debug.print("\t{d} + {d}i, ", .{ z.re, z.im });
        }
        std.debug.print("\n", .{});
    }
    std.debug.print("]\n", .{});

    var E: zml.Array(f64) = try zml.abs(D, .{ .array_allocator = a });
    defer E.deinit(a);

    std.debug.print("E = [", .{});
    for (0..E.shape[0]) |i| {
        std.debug.print("{d}, ", .{E.data[i]});
    }
    std.debug.print("]\n", .{});

    var F = try E.slice(&.{try .init(null, null, -2)});
    std.debug.print("F = E[::-2]\n", .{});
    std.debug.print("F = [", .{});
    for (0..F.shape[0]) |i| {
        std.debug.print("{d}, ", .{(try F.get(&.{i})).*});
    }
    std.debug.print("]\n", .{});

    try zml.ceil_(&F, F, .{});

    std.debug.print("E after ceil = [", .{});
    for (0..E.shape[0]) |i| {
        std.debug.print("{d}, ", .{E.data[i]});
    }
    std.debug.print("]\n", .{});

    std.debug.print("F after ceil = [", .{});
    for (0..F.shape[0]) |i| {
        std.debug.print("{d}, ", .{(try F.get(&.{i})).*});
    }
    std.debug.print("]\n", .{});

    var G: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 5, 5 }, .{});
    defer G.deinit(a);

    for (0..G.shape[0]) |i| {
        for (0..G.shape[1]) |j| {
            try G.set(&.{ i, j }, -zml.scast(f64, (i + 1) * 10 + (j + 1)) / 100);
        }
    }

    std.debug.print("G =\n", .{});
    for (0..G.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..G.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try G.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var H: zml.Array(f64) = try G.slice(&.{ try .init(2, 5, 1), try .init(null, null, -2) });
    defer H.deinit(a);

    std.debug.print("H = G[2:5, ::-2]\n", .{});
    std.debug.print("H =\n", .{});
    for (0..H.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..H.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try H.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    try zml.abs_(&H, H, .{});

    std.debug.print("G after abs =\n", .{});
    for (0..G.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..G.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try G.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    std.debug.print("H after abs =\n", .{});
    for (0..H.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..H.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try H.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    const shapes: []const []const u32 = &.{
        &.{ 15, 3, 1 },
        &.{ 3, 1 },
        &.{ 1, 1, 5 },
    };
    _ = shapes;

    const G_T = try G.transpose(null);
    std.debug.print("G_T = G.transpose()\n", .{});
    std.debug.print("G.shape = [  ", .{});
    for (G.shape[0..G.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});
    std.debug.print("G.strides = [  ", .{});
    for (G.metadata.dense.strides[0..G.ndim]) |stride| {
        std.debug.print("{}  ", .{stride});
    }
    std.debug.print("]\n", .{});
    std.debug.print("G_T.shape = [  ", .{});
    for (G_T.shape[0..G_T.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});
    std.debug.print("G_T.strides = [  ", .{});
    for (G_T.metadata.strided.strides[0..G_T.ndim]) |stride| {
        std.debug.print("{}  ", .{stride});
    }
    std.debug.print("]\n", .{});
    std.debug.print("G_T =\n", .{});
    for (0..G_T.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..G_T.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try G_T.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    const H_T = try H.transpose(null);
    std.debug.print("H_T = H.transpose()\n", .{});
    std.debug.print("H_T =\n", .{});
    for (0..H_T.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..H_T.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try H_T.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }
}

fn addTesting(a: std.mem.Allocator) !void {
    var B: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 5, 1, 8 }, .{ .order = .row_major });
    defer B.deinit(a);
    for (0..B.size) |i| {
        B.data[i] = @floatFromInt(i + 1);
    }
    std.debug.print("B =\n", .{});
    for (0..B.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..B.shape[2]) |j| {
            std.debug.print("{!d:.2}  ", .{(try B.get(&.{ i, 0, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var C: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 1, 8 }, .{ .order = .row_major });
    defer C.deinit(a);
    for (0..C.size) |i| {
        C.data[i] = @floatFromInt(i + 1);
    }
    std.debug.print("\nC =\n", .{});
    for (0..C.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..C.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try C.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var D: zml.Array(f64) = try zml.div(B, C, .{ .array_allocator = a });
    defer D.deinit(a);
    std.debug.print("\nD = B + C\n", .{});
    std.debug.print("D.shape = [  ", .{});
    for (D.shape[0..D.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});
    std.debug.print("\nD =\n", .{});
    for (0..D.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..D.shape[2]) |j| {
            std.debug.print("{!d:.2}  ", .{(try D.get(&.{ i, 0, j })).*});
        }
        std.debug.print("\n", .{});
    }

    const B_T = try (try B.reshape(&.{ 5, 8 })).transpose(null);
    std.debug.print("B_T = B.transpose(null)\n", .{});
    std.debug.print("B_T =\n", .{});
    for (0..B_T.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..B_T.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try B_T.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    const C_T = try C.transpose(null);
    std.debug.print("C_T = C.transpose(null)\n", .{});
    std.debug.print("C_T =\n", .{});
    for (0..C_T.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..C_T.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try C_T.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var E: zml.Array(f64) = try zml.add(B_T, C_T, .{ .array_allocator = a });
    defer E.deinit(a);
    std.debug.print("\nE = B_T + C_T\n", .{});
    std.debug.print("\nE =\n", .{});
    for (0..E.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..E.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try E.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var F: zml.Array(f64) = try zml.add(5, E, .{ .array_allocator = a });
    defer F.deinit(a);
    std.debug.print("\nF = 5 + E\n", .{});
    std.debug.print("\nF =\n", .{});
    for (0..F.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..F.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try F.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var G: zml.Array(f64) = try zml.pow(B, C, .{ .array_allocator = a });
    defer G.deinit(a);
    std.debug.print("\nG = B ** C\n", .{});
    std.debug.print("G =\n", .{});
    for (0..G.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..G.shape[2]) |j| {
            std.debug.print("{!d:.2}  ", .{(try G.get(&.{ i, 0, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var H: zml.Array(f64) = try F.slice(&.{ try .init(2, 5, 1), .all_reverse });
    std.debug.print("\nH = F[2:5, ::-1]\n", .{});
    std.debug.print("H =\n", .{});
    for (0..H.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..H.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try H.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var I: zml.Array(f64) = try .init(a, &.{ H.shape[0], H.shape[1] }, .{});
    defer I.deinit(a);
    for (0..I.size) |i| {
        I.data[i] = @floatFromInt(i + 1);
    }
    std.debug.print("\nI =\n", .{});
    for (0..I.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..I.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try I.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    try zml.ge_(&H, H, I, .{});
    std.debug.print("\nH after add(I) =\n", .{});
    for (0..H.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..H.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try H.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    std.debug.print("\nF after add(I) =\n", .{});
    for (0..F.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..F.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try F.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    const sum: f64 = zml.linalg.blas.dasum(zml.scast(i32, F.size), F.data.ptr, 1);
    std.debug.print("Sum of F: {d}\n", .{sum});
}

fn iterTesting(a: std.mem.Allocator) !void {
    var iterArrR: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 3, 2, 4 }, .{ .order = .row_major });
    defer iterArrR.deinit(a);
    var iterArrC: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 3, 2, 4 }, .{ .order = .col_major });
    defer iterArrC.deinit(a);

    var iterR: zml.array.Iterator(f64) = zml.array.Iterator(f64).init(&iterArrR);
    var iterC: zml.array.Iterator(f64) = zml.array.Iterator(f64).init(&iterArrC);
    std.debug.print("Position(R)    Position(C)     R   C\n", .{});
    std.debug.print("----------------------\n", .{});
    std.debug.print("[  ", .{});
    for (0..iterR.arr.ndim) |i| {
        std.debug.print("{}  ", .{iterR.position[i]});
    }
    std.debug.print("]  [  ", .{});
    for (0..iterC.arr.ndim) |i| {
        std.debug.print("{}  ", .{iterC.position[i]});
    }
    std.debug.print("],  {},  {}\n", .{ iterR.index, iterC.index });
    while (iterR.nextAO(0, .leftToRight) != null and iterC.nextAO(0, .leftToRight) != null) {
        std.debug.print("[  ", .{});
        for (0..iterR.arr.ndim) |i| {
            std.debug.print("{}  ", .{iterR.position[i]});
        }
        std.debug.print("]  [  ", .{});
        for (0..iterC.arr.ndim) |i| {
            std.debug.print("{}  ", .{iterC.position[i]});
        }
        std.debug.print("],  {},  {}\n", .{ iterR.index, iterC.index });
    }
    _ = iterC.nextAO(0, .leftToRight); // Once R has reached the end it is null and the if never reaches C, so C has to run again
    std.debug.print("Final state:\n", .{});
    std.debug.print("[  ", .{});
    for (0..iterR.arr.ndim) |i| {
        std.debug.print("{}  ", .{iterR.position[i]});
    }
    std.debug.print("]  [  ", .{});
    for (0..iterC.arr.ndim) |i| {
        std.debug.print("{}  ", .{iterC.position[i]});
    }
    std.debug.print("],  {},  {}\n\n", .{ iterR.index, iterC.index });
}

fn iterPerfTesting(a: std.mem.Allocator) !void {
    std.debug.print("Column major array, long next, release fast\n", .{});

    var arrBig: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 100, 100, 100, 100 }, .{ .order = .col_major });
    defer arrBig.deinit(a);
    var iterBig: zml.array.Iterator(f64) = zml.array.Iterator(f64).init(&arrBig);

    //
    const n: u32 = 10;
    var start_time = std.time.nanoTimestamp();

    var count: u128 = 0;
    for (0..n) |_| {
        while (iterBig.nextAO(iterBig.array.ndim - 1, .rightToLeft) != null) {
            count += 1;
        }
    }

    var end_time = std.time.nanoTimestamp();
    var duration_ns = end_time - start_time;

    // Convert nanoseconds to seconds as a floating-point number.
    var duration_s: f128 = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Row major iteration took: {d:.9} seconds\n", .{duration_s});

    start_time = std.time.nanoTimestamp();

    for (0..n) |_| {
        while (iterBig.nextAO(0, .leftToRight) != null) {
            count += 1;
        }
    }

    end_time = std.time.nanoTimestamp();
    duration_ns = end_time - start_time;

    // Convert nanoseconds to seconds as a floating-point number.
    duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Column major iteration took: {d:.9} seconds\n", .{duration_s});

    start_time = std.time.nanoTimestamp();

    for (0..n) |_| {
        while (iterBig.next() != null) {
            count += 1;
        }
    }

    end_time = std.time.nanoTimestamp();
    duration_ns = end_time - start_time;

    // Convert nanoseconds to seconds as a floating-point number.
    duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Default iteration took: {d:.9} seconds\n", .{duration_s});
}

fn multiIterTesting(a: std.mem.Allocator) !void {
    var arr1: zml.Array(f64) = try zml.Array(f64).init(a, &[_]u32{ 4, 2, 1 }, .{ .order = .RowMajor });
    defer arr1.deinit();
    var arr2: zml.Array(f64) = try zml.Array(f64).init(a, &[_]u32{ 2, 3 }, .{ .order = .ColumnMajor });
    defer arr2.deinit();

    // Other arrays for broadcasting testing.
    //var arr3: zml.Array(f64) = try zml.Array(f64).init(a, &[_]u32{ 2, 1, 1, 3 }, zml.array.Flags{ .RowMajorContiguous = true, .ColumnMajorContiguous = false });
    var arr3: zml.Array(f64) = try zml.Array(f64).init(a, &[_]u32{ 2, 4, 2, 3 }, .{ .order = .RowMajor });
    defer arr3.deinit();
    var scalar: zml.Array(f64) = try zml.Array(f64).init(a, &[_]u32{}, .{});
    defer scalar.deinit();

    var iter: zml.array.MultiIterator(f64) = try zml.array.MultiIterator(f64).init(&[_]zml.Array(f64){ arr1, arr2, arr3, scalar }, .{ .order = .ColumnMajor });

    std.debug.print("iter.shape = [  ", .{});
    for (0..iter.ndim) |i| {
        std.debug.print("{}  ", .{iter.shape[i]});
    }
    std.debug.print("]\n", .{});
    std.debug.print("arr1.shape = [  ", .{});
    for (0..arr1.ndim) |i| {
        std.debug.print("{}  ", .{arr1.shape[i]});
    }
    std.debug.print("]\n", .{});
    std.debug.print("arr2.shape = [  ", .{});
    for (0..arr2.ndim) |i| {
        std.debug.print("{}  ", .{arr2.shape[i]});
    }
    std.debug.print("]\n", .{});
    std.debug.print("arr3.shape = [  ", .{});
    for (0..arr3.ndim) |i| {
        std.debug.print("{}  ", .{arr3.shape[i]});
    }
    std.debug.print("]\n", .{});

    std.debug.print("F                      1                   2                3                      Scalar\n", .{});
    std.debug.print("-------------------------------------------------------------------------------------------\n", .{});
    std.debug.print("[  ", .{});
    for (0..iter.ndim) |i| {
        std.debug.print("{}  ", .{iter.position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.index});
    for (0..iter.iterators[0].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[0].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[0].index});
    for (0..iter.iterators[1].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[1].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[1].index});
    for (0..iter.iterators[2].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[2].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[2].index});
    for (0..iter.iterators[3].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[3].position[i]});
    }
    std.debug.print("]: {:0>2}\n", .{iter.iterators[3].index});
    while (iter.nextOrder(.RowMajor) != null) {
        std.debug.print("[  ", .{});
        for (0..iter.ndim) |i| {
            std.debug.print("{}  ", .{iter.position[i]});
        }
        std.debug.print("]: {:0>2}   [  ", .{iter.index});
        for (0..iter.iterators[0].ndim) |i| {
            std.debug.print("{}  ", .{iter.iterators[0].position[i]});
        }
        std.debug.print("]: {:0>2}   [  ", .{iter.iterators[0].index});
        for (0..iter.iterators[1].ndim) |i| {
            std.debug.print("{}  ", .{iter.iterators[1].position[i]});
        }
        std.debug.print("]: {:0>2}   [  ", .{iter.iterators[1].index});
        for (0..iter.iterators[2].ndim) |i| {
            std.debug.print("{}  ", .{iter.iterators[2].position[i]});
        }
        std.debug.print("]: {:0>2}   [  ", .{iter.iterators[2].index});
        for (0..iter.iterators[3].ndim) |i| {
            std.debug.print("{}  ", .{iter.iterators[3].position[i]});
        }
        std.debug.print("]: {:0>2}\n", .{iter.iterators[3].index});
    }
    std.debug.print("Final state:\n", .{});
    std.debug.print("[  ", .{});
    for (0..iter.ndim) |i| {
        std.debug.print("{}  ", .{iter.position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.index});
    for (0..iter.iterators[0].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[0].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[0].index});
    for (0..iter.iterators[1].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[1].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[1].index});
    for (0..iter.iterators[2].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[2].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[2].index});
    for (0..iter.iterators[3].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[3].position[i]});
    }
    std.debug.print("]: {:0>2}\n", .{iter.iterators[3].index});
}

fn formatValueWithCustomPadding(value: i32, padding: u8) []const u8 {
    const buffer = std.mem.Allocator.buffer(@sizeOf(value), padding);
    const fmt = "{02d}";
    const len = std.fmt.bufPrint(buffer[0..], fmt, .{value}) catch return "";
    return buffer[0..len];
}
