const std = @import("std");
const zsl = @import("zsl");

pub fn main(init: std.process.Init) !void {
    @setEvalBranchQuota(100_000);

    const N: type = zsl.Complex(zsl.Dyadic(256, 32));
    const layout: zsl.matrix.Layout = .row_major;

    // const arena = init.arena.allocator();
    const gpa = init.gpa;

    const io = init.io;

    var xoshiro = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(io).toMicroseconds()));
    const prng = xoshiro.random();
    const uniform = zsl.stats.Uniform(N).init(zsl.numeric.cast(N, -1), zsl.numeric.cast(N, 1));

    // const m = 6;
    const n = 50;

    var x: zsl.matrix.general.Dense(N, layout) = try .initFn(gpa, n, n, zsl.stats.Uniform(N).sample, .{ uniform, prng });
    defer x.deinit(gpa);

    var xxh: zsl.matrix.general.Dense(N, layout) = try zsl.linalg.matmulAlloc(gpa, x, x.adjointView());
    defer xxh.deinit(gpa);

    const a: zsl.matrix.hermitian.Dense(N, .lower, layout) = try xxh.hermitianView(.lower);

    // std.debug.print("A: {f}\n", .{a.formatter("{d:.4}")});

    var utu: zsl.linalg.utu.Dense(N, layout) = try .init(gpa, n);
    defer utu.deinit(gpa);
    try zsl.linalg.utu.factorInto(&utu, a);

    var a_recreated = try zsl.linalg.matmulAlloc(gpa, utu.u().adjointView(), utu.u());
    defer a_recreated.deinit(gpa);

    var diff = try zsl.matrix.subAlloc(gpa, a, a_recreated);
    defer diff.deinit(gpa);

    std.debug.print("‖A - Uᴴ * U‖ = {e:.4}\n", .{zsl.numeric.cast(f64, try zsl.linalg.normAlloc(gpa, diff, .frobenius))});

    // var a: zsl.matrix.general.Dense(f64, .col_major) = try .initFn(gpa, m, n, zsl.stats.Normal(f64).sample, .{ normal, prng });
    // defer a.deinit(gpa);

    // var a_clone = try a.clone(gpa);
    // defer a_clone.deinit(gpa);

    // var s: zsl.matrix.diagonal.Sparse(f64) = try .init(gpa, m, n);
    // defer s.deinit(gpa);
    // var u: zsl.matrix.general.Dense(f64, .col_major) = try .init(gpa, m, m);
    // defer u.deinit(gpa);
    // var vt: zsl.matrix.general.Dense(f64, .col_major) = try .init(gpa, n, n);
    // defer vt.deinit(gpa);

    // const superb_len = zsl.int.max(1, zsl.int.min(m, n));
    // const superb = try gpa.alloc(f64, superb_len);
    // defer gpa.free(superb);

    // const info = zsl.linalg.lapacke.dgesvd(
    //     zsl.linalg.cblas.layout.col_major,
    //     zsl.linalg.lapacke.job.all,
    //     zsl.linalg.lapacke.job.all,
    //     m,
    //     n,
    //     a_clone.data,
    //     zsl.numeric.cast(isize, a_clone.ld),
    //     s.data,
    //     u.data,
    //     zsl.numeric.cast(isize, u.ld),
    //     vt.data,
    //     zsl.numeric.cast(isize, vt.ld),
    //     superb.ptr,
    // );

    // if (info != 0)
    //     return error.SVD;

    // var us = try zsl.linalg.matmulAlloc(gpa, u, s);
    // defer us.deinit(gpa);
    // var a_reconstructed = try zsl.linalg.matmulAlloc(gpa, us, vt);
    // defer a_reconstructed.deinit(gpa);

    // printMatrix("A = U S V^T", a);
    // printMatrix("U", u);
    // printMatrix("S", s);
    // printMatrix("V^T", vt);
    // printMatrix("A = U S V^T (reconstructed)", a_reconstructed);

    // var diff = try zsl.matrix.subAlloc(gpa, a_reconstructed, a);
    // defer diff.deinit(gpa);

    // const diff_norm = try zsl.linalg.normAlloc(gpa, diff, .frobenius);

    // std.debug.print("‖A - U S V^T‖ = {e}\n", .{diff_norm});
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

/// Formats the dyadic exactly into a base-10 string.
/// The caller takes ownership of the returned slice and must free it.
pub fn dyadicToString(allocator: std.mem.Allocator, x: anytype) ![]u8 {
    // Handle edge cases
    if (x.isNan())
        return allocator.dupe(u8, "NaN");

    if (x.isInf())
        return allocator.dupe(u8, if (x.positive) "Inf" else "-Inf");

    if (x.isZero())
        return allocator.dupe(u8, if (x.positive) "0.0" else "-0.0");

    const Mantissa = @TypeOf(x).Mantissa;
    const bits = @typeInfo(Mantissa).int.bits;

    // Calculate max significant figures at compile time
    const calculated_sig_figs = comptime @as(usize, @intFromFloat(@as(f64, @floatFromInt(bits)) * 0.3010299956639812 * 1.05));
    const max_sig_figs: usize = if (calculated_sig_figs == 0) 1 else calculated_sig_figs;

    const Helpers = struct {
        // Add two base-10 digit arrays
        fn addDigits(alloc: std.mem.Allocator, dest: *std.ArrayList(u8), src: []const u8) !void {
            var carry: u8 = 0;
            var i: usize = 0;
            while (i < src.len or carry > 0) : (i += 1) {
                if (i >= dest.items.len) {
                    try dest.append(alloc, 0);
                }
                const a = dest.items[i];
                const b = if (i < src.len) src[i] else 0;
                const sum = a + b + carry;
                dest.items[i] = sum % 10;
                carry = sum / 10;
            }
        }

        // Multiply a base-10 digit array by a small scalar
        fn mulByDigit(alloc: std.mem.Allocator, dest: *std.ArrayList(u8), multiplier: u32) !void {
            var carry: u32 = 0;
            for (dest.items) |*digit| {
                const prod = @as(u32, digit.*) * multiplier + carry;
                digit.* = @intCast(prod % 10);
                carry = prod / 10;
            }
            while (carry > 0) {
                try dest.append(alloc, @intCast(carry % 10));
                carry /= 10;
            }
        }
    };

    var num: std.ArrayList(u8) = .empty;
    defer num.deinit(allocator);
    try num.append(allocator, 0);

    var power_of_2: std.ArrayList(u8) = .empty;
    defer power_of_2.deinit(allocator);
    try power_of_2.append(allocator, 1);

    // Step 1: Load the base mantissa into the `num` bignum array
    var i: u16 = 0;
    while (i < bits) : (i += 1) {
        const shift: std.math.Log2Int(Mantissa) = @intCast(i);
        const bit = (x.mantissa >> shift) & 1;
        if (bit == 1) {
            try Helpers.addDigits(allocator, &num, power_of_2.items);
        }
        try Helpers.mulByDigit(allocator, &power_of_2, 2);
    }

    // Step 2: Apply the exponent
    const e_i32: i32 = x.exponent;
    const abs_e: u32 = @intCast(if (e_i32 < 0) -e_i32 else e_i32);

    if (e_i32 > 0) {
        // Number is Mantissa * 2^E
        var j: u32 = 0;
        while (j < abs_e) : (j += 1) {
            try Helpers.mulByDigit(allocator, &num, 2);
        }
    } else if (e_i32 < 0) {
        // Number is Mantissa * 5^|E| / 10^|E|
        // We compute the numerator exactly to avoid floating-point logic
        var j: u32 = 0;
        while (j < abs_e) : (j += 1) {
            try Helpers.mulByDigit(allocator, &num, 5);
        }
    }

    // Pad with trailing zeros to guarantee we safely cross the decimal point
    if (e_i32 < 0) {
        while (num.items.len <= abs_e) {
            try num.append(allocator, 0);
        }
    }

    // Step 3: Format the final string
    var result: std.ArrayList(u8) = .empty;
    if (!x.positive) {
        try result.append(allocator, '-');
    }

    var sig_figs: usize = 0;
    var started_counting = false;

    if (e_i32 >= 0) {
        // Integer representation
        var idx: usize = num.items.len;
        while (idx > 0) {
            idx -= 1;
            const digit = num.items[idx];

            if (digit != 0) started_counting = true;
            if (started_counting) sig_figs += 1;

            if (sig_figs <= max_sig_figs) {
                try result.append(allocator, digit + '0');
            } else {
                // Cap exceeded: Pad integer part with zeros for magnitude safety
                try result.append(allocator, '0');
            }
        }
        try result.appendSlice(allocator, ".0");
    } else {
        // Fractional representation, inject decimal point
        var idx: usize = num.items.len;
        while (idx > 0) {
            idx -= 1;
            const digit = num.items[idx];

            if (digit != 0 and !started_counting) started_counting = true;
            if (started_counting) sig_figs += 1;

            if (sig_figs <= max_sig_figs or !started_counting) {
                try result.append(allocator, digit + '0');
            } else {
                // Cap exceeded
                if (idx >= abs_e) {
                    // Haven't hit the decimal yet, pad to preserve scale
                    try result.append(allocator, '0');
                } else {
                    // Safely right of the decimal point, just truncate
                    break;
                }
            }

            if (idx == abs_e) {
                try result.append(allocator, '.');
            }
        }

        // Cleanup: If the cap caused the string to stop exactly on the decimal, append a '0'
        if (result.items.len > 0 and result.items[result.items.len - 1] == '.') {
            try result.append(allocator, '0');
        }
    }

    return result.toOwnedSlice(allocator);
}

// /// Generate a random m×n matrix A with specified 2-norm condition number `kappa`.
// pub fn random_matrix(
//     allocator: std.mem.Allocator,
//     m: u32,
//     n: u32,
//     kappa: f64,
//     order: zsl.Layout,
// ) ![]f64 {
//     // - `allocator`: memory allocator
//     // - `m`, `n`: dimensions
//     // - `kappa`: desired condition number κ₂(A)
//     // - `order`: either .RowMajor or .ColMajor
//     //
//     // Method:
//     // 1. Construct singular values σᵢ geometrically spaced between 1 and κ^(1/min(m,n)).
//     // 2. Generate random m×min(m,n) matrix X, QR factor → Q₁ (m×r).
//     // 3. Generate random n×min(m,n) matrix Y, QR factor → Q₂ (n×r).
//     // 4. Form A = Q₁ * diag(σ) * Q₂ᵀ.
//     //     - Compute B = diag(σ) * Q₂ᵀ.
//     //     - Then A = Q₁ × B.
//     const r = zsl.int.min(m, n);

//     // 1) geometric singular values σ[i]
//     var sig = try allocator.alloc(f64, r);
//     defer allocator.free(sig);
//     for (0..r) |i| {
//         sig[i] = zsl.float.pow(kappa, zsl.float.div(zsl.scast(f64, i), r - 1));
//     }

//     // 2) random X: m×r → QR → Q₁
//     const X = try random_buffer(allocator, m * r);
//     defer allocator.free(X);
//     const tauX = try allocator.alloc(f64, r);
//     defer allocator.free(tauX);
//     _ = ci.LAPACKE_dgeqrf(
//         order.toCInt(),
//         zsl.scast(c_int, m),
//         zsl.scast(c_int, r),
//         X.ptr,
//         zsl.scast(c_int, if (order == .col_major) m else r),
//         tauX.ptr,
//     );
//     _ = ci.LAPACKE_dorgqr(
//         order.toCInt(),
//         zsl.scast(c_int, m),
//         zsl.scast(c_int, r),
//         zsl.scast(c_int, r),
//         X.ptr,
//         zsl.scast(c_int, if (order == .col_major) m else r),
//         tauX.ptr,
//     );
//     // X now holds Q₁

//     // 3) random Y: n×r → QR → Q₂
//     const Y = try random_buffer(allocator, n * r);
//     defer allocator.free(Y);
//     const tauY = try allocator.alloc(f64, r);
//     defer allocator.free(tauY);
//     _ = ci.LAPACKE_dgeqrf(
//         order.toCInt(),
//         zsl.scast(c_int, n),
//         zsl.scast(c_int, r),
//         Y.ptr,
//         zsl.scast(c_int, if (order == .col_major) n else r),
//         tauY.ptr,
//     );
//     _ = ci.LAPACKE_dorgqr(
//         order.toCInt(),
//         zsl.scast(c_int, n),
//         zsl.scast(c_int, r),
//         zsl.scast(c_int, r),
//         Y.ptr,
//         zsl.scast(c_int, if (order == .col_major) n else r),
//         tauY.ptr,
//     );
//     // Y now holds Q₂

//     // 4a) form B = diag(sig) * Q2ᵀ (B is r×n in same layout)
//     const B = try allocator.alloc(f64, r * n);
//     defer allocator.free(B);
//     for (0..r) |i| {
//         const row_start = if (order == .col_major)
//             B.ptr + i
//         else
//             B.ptr + i * n;

//         const stride = if (order == .col_major) r else 1;

//         for (0..n) |j| {
//             // Q2[j,i] is at Y.ptr + j*ldY + i
//             const q2ji = (Y.ptr + j * if (order == .col_major) n else r)[i];
//             row_start[j * stride] = sig[i] * q2ji;
//         }
//     }

//     // 4b) A = Q1 (m×r) × B (r×n) → A (m×n)
//     const A = try allocator.alloc(f64, m * n);
//     zsl.linalg.blas.dgemm(
//         order,
//         .no_trans,
//         .no_trans,
//         zsl.scast(i32, m),
//         zsl.scast(i32, n),
//         zsl.scast(i32, r),
//         1,
//         X.ptr,
//         zsl.scast(i32, if (order == .col_major) m else r),
//         B.ptr,
//         zsl.scast(i32, if (order == .col_major) r else n),
//         0,
//         A.ptr,
//         zsl.scast(i32, if (order == .col_major) m else n),
//     );

//     return A;
// }

// pub fn random_matrix_buffer(
//     allocator: std.mem.Allocator,
//     m: u32,
//     n: u32,
//     kappa: f64,
//     A: []f64,
//     order: zsl.Order,
// ) !void {
//     // - `allocator`: memory allocator
//     // - `m`, `n`: dimensions
//     // - `kappa`: desired condition number κ₂(A)
//     // - `order`: either .RowMajor or .ColMajor
//     //
//     // Method:
//     // 1. Construct singular values σᵢ geometrically spaced between 1 and κ^(1/min(m,n)).
//     // 2. Generate random m×min(m,n) matrix X, QR factor → Q₁ (m×r).
//     // 3. Generate random n×min(m,n) matrix Y, QR factor → Q₂ (n×r).
//     // 4. Form A = Q₁ * diag(σ) * Q₂ᵀ.
//     //     - Compute B = diag(σ) * Q₂ᵀ.
//     //     - Then A = Q₁ × B.
//     const r = zsl.int.min(m, n);

//     // 1) geometric singular values σ[i]
//     var sig = try allocator.alloc(f64, r);
//     defer allocator.free(sig);
//     for (0..r) |i| {
//         sig[i] = zsl.float.pow(kappa, zsl.float.div(zsl.scast(f64, i), r - 1));
//     }

//     // 2) random X: m×r → QR → Q₁
//     const X = try random_buffer(allocator, m * r);
//     defer allocator.free(X);
//     const tauX = try allocator.alloc(f64, r);
//     defer allocator.free(tauX);
//     _ = ci.LAPACKE_dgeqrf(
//         order.toCInt(),
//         zsl.scast(c_int, m),
//         zsl.scast(c_int, r),
//         X.ptr,
//         zsl.scast(c_int, if (order == .col_major) m else r),
//         tauX.ptr,
//     );
//     _ = ci.LAPACKE_dorgqr(
//         order.toCInt(),
//         zsl.scast(c_int, m),
//         zsl.scast(c_int, r),
//         zsl.scast(c_int, r),
//         X.ptr,
//         zsl.scast(c_int, if (order == .col_major) m else r),
//         tauX.ptr,
//     );
//     // X now holds Q₁

//     // 3) random Y: n×r → QR → Q₂
//     const Y = try random_buffer(allocator, n * r);
//     defer allocator.free(Y);
//     const tauY = try allocator.alloc(f64, r);
//     defer allocator.free(tauY);
//     _ = ci.LAPACKE_dgeqrf(
//         order.toCInt(),
//         zsl.scast(c_int, n),
//         zsl.scast(c_int, r),
//         Y.ptr,
//         zsl.scast(c_int, if (order == .col_major) n else r),
//         tauY.ptr,
//     );
//     _ = ci.LAPACKE_dorgqr(
//         order.toCInt(),
//         zsl.scast(c_int, n),
//         zsl.scast(c_int, r),
//         zsl.scast(c_int, r),
//         Y.ptr,
//         zsl.scast(c_int, if (order == .col_major) n else r),
//         tauY.ptr,
//     );
//     // Y now holds Q₂

//     // 4a) form B = diag(sig) * Q2ᵀ (B is r×n in same layout)
//     const B = try allocator.alloc(f64, r * n);
//     defer allocator.free(B);
//     for (0..r) |i| {
//         const row_start = if (order == .col_major)
//             B.ptr + i
//         else
//             B.ptr + i * n;

//         const stride = if (order == .col_major) r else 1;

//         for (0..n) |j| {
//             // Q2[j,i] is at Y.ptr + j*ldY + i
//             const q2ji = (Y.ptr + j * if (order == .col_major) n else r)[i];
//             row_start[j * stride] = sig[i] * q2ji;
//         }
//     }

//     // 4b) A = Q1 (m×r) × B (r×n) → A (m×n)
//     zsl.linalg.blas.dgemm(
//         order,
//         .no_trans,
//         .no_trans,
//         zsl.scast(i32, m),
//         zsl.scast(i32, n),
//         zsl.scast(i32, r),
//         1,
//         X.ptr,
//         zsl.scast(i32, if (order == .col_major) m else r),
//         B.ptr,
//         zsl.scast(i32, if (order == .col_major) r else n),
//         0,
//         A.ptr,
//         zsl.scast(i32, if (order == .col_major) m else n),
//     );
// }

// /// Generates a random symmetric positive definite matrix with a specified condition number.
// fn random_symmetric_positive_definite_matrix(
//     allocator: std.mem.Allocator,
//     size: u32,
//     cond_target: f64,
// ) ![]f64 {
//     // 1) compute λ_i geometrically between 1 and cond_target
//     var lambdas = try allocator.alloc(f64, size);
//     defer allocator.free(lambdas);
//     const l_min = 1;
//     const l_max = cond_target;
//     for (0..size) |i| {
//         lambdas[i] = l_min * zsl.float.pow(l_max / l_min, zsl.float.div(zsl.scast(f64, i), size - 1));
//     }

//     // 2) random X, then QR -> get Q (n×n)
//     const X = try random_buffer(allocator, size * size);
//     defer allocator.free(X);
//     const tau = try allocator.alloc(f64, size);
//     defer allocator.free(tau);
//     // QR factorization in-place: X = Q*R
//     _ = ci.LAPACKE_dgeqrf(
//         zsl.Order.col_major.toCInt(),
//         zsl.scast(c_int, size),
//         zsl.scast(c_int, size),
//         X.ptr,
//         zsl.scast(c_int, size),
//         tau.ptr,
//     );
//     _ = ci.LAPACKE_dorgqr(
//         zsl.Order.col_major.toCInt(),
//         zsl.scast(c_int, size),
//         zsl.scast(c_int, size),
//         zsl.scast(c_int, size),
//         X.ptr,
//         zsl.scast(c_int, size),
//         tau.ptr,
//     );
//     // now X holds Q

//     // 3) form B = D * Q^T, where D = diag(λ_1, λ_2, ..., λ_n)
//     var B = try allocator.alloc(f64, size * size);
//     defer allocator.free(B);
//     // copy Q^T into B
//     for (0..size) |i| {
//         for (0..size) |j|
//             B[i * size + j] = X[j * size + i];
//     }
//     // scale row i of B by λ_i
//     for (0..size) |i| {
//         zsl.linalg.blas.dscal(zsl.scast(i32, size), lambdas[i], B.ptr + i * size, 1);
//     }

//     // 4) Assemble A = Q * D * Q^T as A = Q * B
//     const A = try allocator.alloc(f64, size * size);
//     zsl.linalg.blas.dgemm(
//         .col_major,
//         .no_trans,
//         .no_trans,
//         zsl.scast(i32, size),
//         zsl.scast(i32, size),
//         zsl.scast(i32, size),
//         1,
//         X.ptr,
//         zsl.scast(i32, size),
//         B.ptr,
//         zsl.scast(i32, size),
//         0,
//         A.ptr,
//         zsl.scast(i32, size),
//     );

//     return A;
// }

// fn frobernius_norm_difference_matrix(a: anytype, b: anytype) !f64 {
//     const m = if (comptime zsl.meta.isSymmetricMatrix(@TypeOf(a)) or
//         zsl.meta.isHermitianMatrix(@TypeOf(a)) or
//         zsl.meta.isTridiagonalMatrix(@TypeOf(a)) or
//         zsl.meta.isPermutationMatrix(@TypeOf(a)))
//         a.size
//     else
//         a.rows;
//     const n = if (comptime zsl.meta.isSymmetricMatrix(@TypeOf(a)) or
//         zsl.meta.isHermitianMatrix(@TypeOf(a)) or
//         zsl.meta.isTridiagonalMatrix(@TypeOf(a)) or
//         zsl.meta.isPermutationMatrix(@TypeOf(a)))
//         a.size
//     else
//         a.cols;

//     var norm: f64 = 0;

//     var i: u32 = 0;
//     while (i < m) : (i += 1) {
//         var j: u32 = 0;
//         while (j < n) : (j += 1) {
//             if (comptime !zsl.meta.isComplex(zsl.meta.Numeric(@TypeOf(a))) and !zsl.meta.isComplex(zsl.meta.Numeric(@TypeOf(b)))) {
//                 const diff = try a.get(i, j) - try b.get(i, j);
//                 norm += diff * diff;
//             } else {
//                 const diff = try zsl.sub(try a.get(i, j), try b.get(i, j), .{});
//                 norm += try zsl.abs2(diff, .{});
//             }
//         }
//     }

//     return zsl.float.sqrt(norm);
// }

// fn random_complex_hermitian_positive_definite_matrix(
//     allocator: std.mem.Allocator,
//     size: u32,
//     factor: f64,
// ) ![]zsl.cf64 {
//     // allocate M
//     const M = try random_complex_matrix(allocator, size, size, factor);
//     // allocate A = M M^H
//     const A = try allocator.alloc(zsl.cf64, size * size);

//     // A = M × M^H
//     zsl.linalg.blas.zgemm(
//         .row_major,
//         .no_trans,
//         .conj_trans,
//         zsl.scast(i32, size),
//         zsl.scast(i32, size),
//         zsl.scast(i32, size),
//         zsl.cf64.init(1, 0),
//         M.ptr,
//         zsl.scast(i32, size),
//         M.ptr,
//         zsl.scast(i32, size),
//         zsl.cf64.init(0, 0),
//         A.ptr,
//         zsl.scast(i32, size),
//     );

//     allocator.free(M);
//     return A;
// }
