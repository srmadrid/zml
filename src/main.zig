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

fn random_buffer_fill(
    buffer: []f64,
) void {
    var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
    //var prng = std.Random.DefaultPrng.init(2); // fixed seed for reproducibility
    const rand = prng.random();

    for (0..buffer.len) |i| {
        buffer[i] = rand.float(f64);
    }
}

fn random_buffer_fill_complex(
    buffer: []zml.cf64,
) void {
    //var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
    var prng = std.Random.DefaultPrng.init(2); // fixed seed for reproducibility
    const rand = prng.random();

    for (0..buffer.len) |i| {
        buffer[i] = zml.cf64.init(rand.float(f64), rand.float(f64));
    }
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

pub fn random_matrix_buffer(
    allocator: std.mem.Allocator,
    m: u32,
    n: u32,
    kappa: f64,
    A: []f64,
    order: zml.Order,
) !void {
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

fn frobernius_norm_difference_matrix(a: anytype, b: anytype) !f64 {
    const m = if (comptime zml.types.isSymmetricMatrix(@TypeOf(a)) or
        zml.types.isHermitianMatrix(@TypeOf(a)) or
        zml.types.isTridiagonalMatrix(@TypeOf(a)) or
        zml.types.isPermutationMatrix(@TypeOf(a)))
        a.size
    else
        a.rows;
    const n = if (comptime zml.types.isSymmetricMatrix(@TypeOf(a)) or
        zml.types.isHermitianMatrix(@TypeOf(a)) or
        zml.types.isTridiagonalMatrix(@TypeOf(a)) or
        zml.types.isPermutationMatrix(@TypeOf(a)))
        a.size
    else
        a.cols;

    var norm: f64 = 0;

    var i: u32 = 0;
    while (i < m) : (i += 1) {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            if (comptime !zml.types.isComplex(zml.types.Numeric(@TypeOf(a))) and !zml.types.isComplex(zml.types.Numeric(@TypeOf(b)))) {
                const diff = try a.get(i, j) - try b.get(i, j);
                norm += diff * diff;
            } else {
                const diff = try zml.sub(try a.get(i, j), try b.get(i, j), .{});
                norm += try zml.abs2(diff, .{});
            }
        }
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

    // try binopPerfTesting(a);

    // try decompPerfTesting(a);

    try matrixTesting(a);
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

fn fill(a: anytype, factor: u32) void {
    const A: type = zml.types.Numeric(@TypeOf(a));
    if (comptime zml.types.isVector(@TypeOf(a))) {
        var i: u32 = 0;
        while (i < a.len) : (i += 1) {
            if (comptime zml.types.isComplex(A)) {
                a.data[i] = A.init(zml.scast(zml.types.Scalar(A), i + factor), zml.scast(zml.types.Scalar(A), i + factor));
            } else {
                a.data[i] = zml.scast(A, i + factor);
            }
        }

        return;
    }

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
            while (i < zml.int.min(a.rows, a.cols)) : (i += 1) {
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
        .permutation => {
            randomPermutation(a.data[0..a.size]);
        },
        else => unreachable,
    }
}

fn print(name: []const u8, a: anytype) void {
    if (comptime zml.types.isVector(@TypeOf(a))) {
        std.debug.print("Vector {s}:\n", .{name});
        var i: u32 = 0;
        while (i < a.len) : (i += 1) {
            if (comptime zml.types.isComplex(zml.types.Numeric(@TypeOf(a)))) {
                std.debug.print("{d:3} + {d:3}i\n", .{ (a.get(i) catch unreachable).re, (a.get(i) catch unreachable).im });
            } else {
                std.debug.print("{d:3}\n", .{a.get(i) catch unreachable});
            }
        }
    } else {
        std.debug.print("Matrix {s}:\n", .{name});
        if (comptime zml.types.isSymmetricDenseMatrix(@TypeOf(a)) or zml.types.isHermitianDenseMatrix(@TypeOf(a)) or zml.types.isTridiagonalDenseMatrix(@TypeOf(a)) or zml.types.isSymmetricSparseMatrix(@TypeOf(a)) or zml.types.isHermitianSparseMatrix(@TypeOf(a)) or zml.types.isSymmetricBlockSparseMatrix(@TypeOf(a)) or zml.types.isHermitianBlockSparseMatrix(@TypeOf(a))) {
            var i: u32 = 0;
            while (i < a.size) : (i += 1) {
                var j: u32 = 0;
                while (j < a.size) : (j += 1) {
                    if (comptime zml.types.isComplex(zml.types.Numeric(@TypeOf(a)))) {
                        std.debug.print("{d:8.4} + {d:8.4}i  ", .{ (a.get(i, j) catch unreachable).re, (a.get(i, j) catch unreachable).im });
                    } else {
                        std.debug.print("{d:8.4}  ", .{a.get(i, j) catch unreachable});
                    }
                }
                std.debug.print("\n", .{});
            }
            std.debug.print("\n", .{});
        } else if (comptime zml.types.isPermutationSparseMatrix(@TypeOf(a))) {
            var i: u32 = 0;
            while (i < a.size) : (i += 1) {
                var j: u32 = 0;
                while (j < a.size) : (j += 1) {
                    if (comptime zml.types.isComplex(zml.types.Numeric(@TypeOf(a)))) {
                        std.debug.print("{d}  ", .{(a.get(i, j) catch unreachable).re});
                    } else {
                        std.debug.print("{d}  ", .{a.get(i, j) catch unreachable});
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
                        std.debug.print("{d:8.4} + {d:8.4}i  ", .{ (a.get(i, j) catch unreachable).re, (a.get(i, j) catch unreachable).im });
                    } else {
                        std.debug.print("{d:8.4}  ", .{a.get(i, j) catch unreachable});
                    }
                }
                std.debug.print("\n", .{});
            }
            std.debug.print("\n", .{});
        }
    }
}

fn randomPermutation(data: []u32) void {
    //std.Thread.sleep(1000000000);

    var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
    const rand = prng.random();

    // Initialize with identity permutation
    var i: u32 = 0;
    while (i < data.len) : (i += 1) {
        data[i] = i;
    }

    // Shuffle using Fisher-Yates algorithm
    i = zml.scast(u32, data.len - 1);
    while (i > 0) : (i -= 1) {
        const j = rand.intRangeAtMost(u32, 0, i);
        const temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }
}

fn binopPerfTesting(a: std.mem.Allocator) !void {
    const print_mats: bool = true;

    var A: zml.matrix.General(f64, .col_major) = try .init(a, 8, 8);
    defer A.deinit(a);

    fill(A, 1);
    if (print_mats) print("A", A);

    var B: zml.matrix.Permutation(f64) = try .init(a, 8);
    defer B.deinit(a);

    fill(B, 2);
    if (print_mats) print("B", B.transpose());

    const start_time = std.time.nanoTimestamp();
    var C = try zml.mul(A, B.transpose(), .{ .matrix_allocator = a });
    const end_time: i128 = std.time.nanoTimestamp();
    defer C.deinit(a);

    std.debug.print("Took: {d} seconds\n\n", .{zml.float.div(end_time - start_time, 1e9)});

    if (print_mats) print("C = A * B", C);
}

fn decompPerfTesting(a: std.mem.Allocator) !void {
    const print_mats: bool = false;

    var A: zml.matrix.General(zml.cf64, .row_major) = try .init(a, 2000, 1000);
    defer A.deinit(a);

    //fill(A, 2);
    random_buffer_fill_complex(A.data[0 .. A.rows * A.cols]);
    if (print_mats) print("A", A);

    const tau: []zml.cf64 = try a.alloc(zml.cf64, zml.int.min(A.rows, A.cols));
    defer a.free(tau);
    const work: []zml.cf64 = try a.alloc(zml.cf64, A.cols * 32);
    defer a.free(work);

    var start_time = std.time.nanoTimestamp();
    var qr = try zml.linalg.qr(a, A, .{});
    var end_time: i128 = std.time.nanoTimestamp();
    defer qr.deinit(a);

    std.debug.print("Decomposition took: {d} seconds\n\n", .{zml.float.div(end_time - start_time, 1e9)});

    var q = try qr.q(a, .full, .{});
    defer q.deinit(a);
    const r = qr.r(.full);

    if (print_mats) print("Q", q);
    if (print_mats) print("R", r);

    start_time = std.time.nanoTimestamp();
    var A_reconstructed = try zml.mul(q, r, .{ .matrix_allocator = a });
    defer A_reconstructed.deinit(a);
    end_time = std.time.nanoTimestamp();

    std.debug.print("Reconstruction took: {d} seconds\n\n", .{zml.float.div(end_time - start_time, 1e9)});

    if (print_mats) print("A reconstructed", A_reconstructed);

    std.debug.print("||A - A_r||_F = {d}\n", .{try frobernius_norm_difference_matrix(A, A_reconstructed)});
}

fn is_square(comptime T: type) bool {
    return zml.types.isSymmetricDenseMatrix(T) or
        zml.types.isHermitianDenseMatrix(T) or
        zml.types.isTridiagonalDenseMatrix(T) or
        zml.types.isSymmetricSparseMatrix(T) or
        zml.types.isHermitianSparseMatrix(T) or
        zml.types.isPermutationSparseMatrix(T);
}

fn random_matrix_t(
    comptime T: type,
    allocator: std.mem.Allocator,
    rows: u32,
    cols: u32,
) !T {
    var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
    const rand = prng.random();

    if (comptime zml.types.isDenseMatrix(T)) {
        var result: T = if (comptime is_square(T))
            try T.init(allocator, rows)
        else
            try T.init(allocator, rows, cols);
        errdefer result.deinit(allocator);

        var i: u32 = 0;
        while (i < rows) : (i += 1) {
            var j: u32 = 0;
            while (j < cols) : (j += 1) {
                if (comptime zml.types.isComplex(zml.types.Numeric(T))) {
                    if (zml.types.isHermitianDenseMatrix(T) and i == j) {
                        // diagonal of Hermitian matrix must be real
                        result.set(i, j, zml.types.Numeric(T).init(rand.float(f64), 0)) catch unreachable;

                        continue;
                    }

                    result.set(i, j, zml.types.Numeric(T).init(rand.float(f64), rand.float(f64))) catch unreachable;
                } else {
                    result.set(i, j, rand.float(zml.types.Numeric(T))) catch unreachable;
                }
            }
        }

        return result;
    } else if (comptime zml.types.isPermutationSparseMatrix(T)) {
        return error.UnsupportedMatrixType;
    } else {}
}

fn print_matrix(desc: []const u8, A: anytype) void {
    std.debug.print("\nMatrix {s}:\n\n", .{desc});

    const rows = if (comptime is_square(@TypeOf(A)))
        A.size
    else
        A.rows;
    const cols = if (comptime is_square(@TypeOf(A)))
        A.size
    else
        A.cols;

    var i: u32 = 0;
    while (i < rows) : (i += 1) {
        std.debug.print("\t", .{});

        var j: u32 = 0;
        while (j < cols) : (j += 1) {
            if (comptime zml.types.isComplex(zml.types.Numeric(f64))) {
                std.debug.print("{d:7.4} + {d:7.4}i  ", .{ (A.get(i, j) catch unreachable).re, (A.get(i, j) catch unreachable).im });
            } else {
                std.debug.print("{d:5.4}  ", .{A.get(i, j) catch unreachable});
            }
        }
        std.debug.print("\n", .{});
    }
    std.debug.print("\n", .{});
}

fn matrixTesting(a: std.mem.Allocator) !void {
    // var A = try random_matrix_t(
    //     zml.matrix.dense.Tridiagonal(f64),
    //     a,
    //     5,
    //     5,
    // );
    var A: zml.matrix.sparse.Builder(f64, .col_major) = try .init(a, 5, 6, 10);
    //defer A.deinit(a);

    try A.set(a, 0, 0, 1);
    try A.set(a, 1, 0, 2);
    try A.set(a, 4, 0, 3);
    try A.set(a, 0, 1, 1);
    try A.set(a, 3, 2, 4);
    try A.set(a, 1, 3, 5);
    try A.set(a, 2, 3, 6);
    try A.set(a, 2, 1, 8); // Try self-ordering
    try A.set(a, 4, 4, 6);
    try A.set(a, 2, 5, 7);

    print_matrix("A", A);

    std.debug.print("A.data = {any}\n", .{A.data[0..A.nnz]});
    std.debug.print("A.row = {any}\n", .{A.row[0..A.nnz]});
    std.debug.print("A.col = {any}\n", .{A.col[0..A.nnz]});

    var B = try A.compile(a);
    defer B.deinit(a);

    print_matrix("B", B);

    std.debug.print("B.data = {any}\n", .{B.data[0..B.nnz]});
    std.debug.print("B.idx = {any}\n", .{B.idx[0..B.nnz]});
    std.debug.print("B.ptr = {any}\n", .{B.ptr[0 .. B.cols + 1]});
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
