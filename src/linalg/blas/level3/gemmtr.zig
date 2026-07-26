const std = @import("std");

const options = @import("options");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

pub fn gemmtr(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    transa: linalg.blas.Transpose,
    transb: linalg.blas.Transpose,
    n: usize,
    k: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    b: anytype,
    ldb: usize,
    beta: anytype,
    c: anytype,
    ldc: usize,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);
    const Be: type = @TypeOf(beta);
    comptime var C: type = @TypeOf(c);

    comptime if (!meta.isNumeric(Al) or !meta.isNumeric(Be) or
        !meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(B) or !meta.isNumeric(meta.Child(B)) or
        !meta.isManyItemPointer(C) or meta.isConstPointer(C) or !meta.isNumeric(meta.Child(C)))
        @compileError("zsl.linalg.blas.gemmtr: alpha and beta must be numerics, a and b must be many-item pointers to numerics, and c must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tb: " ++ @typeName(B) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\tc: " ++ @typeName(C) ++ "\n");

    A = meta.Child(A);
    B = meta.Child(B);
    C = meta.Child(C);

    const notransa = transa == .no_trans or transa == .conj_no_trans;
    const notransb = transb == .no_trans or transb == .conj_no_trans;
    const noconja = transa == .no_trans or transa == .trans;
    const noconjb = transb == .no_trans or transb == .trans;

    if (layout == .col_major) {
        const nrowa: usize = if (notransa) n else k;
        const nrowb: usize = if (notransb) k else n;

        if (lda < int.max(1, nrowa) or ldb < int.max(1, nrowb) or ldc < int.max(1, n))
            return linalg.blas.Error.InvalidArgument;

        // Quick return if possible.
        if (n == 0)
            return;

        if (numeric.eq(alpha, 0)) {
            if (numeric.ne(beta, 1)) {
                for (0..n) |j| {
                    linalg.blas.scal(
                        if (uplo == .upper)
                            j + 1
                        else
                            n - j,
                        beta,
                        c + if (uplo == .upper)
                            j * ldc
                        else
                            j + j * ldc,
                        1,
                    ) catch unreachable;
                }
            }

            return;
        }

        switch (noconja) {
            true => switch (noconjb) {
                true => return k_gemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc, true, true),
                false => return k_gemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc, true, false),
            },
            false => switch (noconjb) {
                true => return k_gemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc, false, true),
                false => return k_gemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc, false, false),
            },
        }
    } else {
        const nrowa: usize = if (notransa) k else n;
        const nrowb: usize = if (notransb) n else k;

        if (lda < int.max(1, nrowa) or ldb < int.max(1, nrowb) or ldc < int.max(1, n))
            return linalg.blas.Error.InvalidArgument;

        // Quick return if possible.
        if (n == 0)
            return;

        if (numeric.eq(alpha, 0)) {
            if (numeric.ne(beta, 1)) {
                for (0..n) |i| {
                    linalg.blas.scal(
                        if (uplo == .upper)
                            n - i
                        else
                            i + 1,
                        beta,
                        c +
                            if (uplo == .upper)
                                i + i * ldc
                            else
                                i * ldc,
                        1,
                    ) catch unreachable;
                }
            }

            return;
        }

        switch (noconja) {
            true => switch (noconjb) {
                true => return k_gemmtr(uplo.invert(), transb, transa, n, k, alpha, b, ldb, a, lda, beta, c, ldc, true, true),
                false => return k_gemmtr(uplo.invert(), transb, transa, n, k, alpha, b, ldb, a, lda, beta, c, ldc, false, true),
            },
            false => switch (noconjb) {
                true => return k_gemmtr(uplo.invert(), transb, transa, n, k, alpha, b, ldb, a, lda, beta, c, ldc, true, false),
                false => return k_gemmtr(uplo.invert(), transb, transa, n, k, alpha, b, ldb, a, lda, beta, c, ldc, false, false),
            },
        }
    }
}

fn k_gemmtr(
    uplo: matrix.Uplo,
    transa: linalg.blas.Transpose,
    transb: linalg.blas.Transpose,
    n: usize,
    k: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    b: anytype,
    ldb: usize,
    beta: anytype,
    c: anytype,
    ldc: usize,
    comptime noconja: bool,
    comptime noconjb: bool,
) !void {
    const A: type = meta.Child(@TypeOf(a));
    const B: type = meta.Child(@TypeOf(b));

    // First form  C = beta * C.
    if (numeric.ne(beta, 1)) {
        for (0..n) |j| {
            linalg.blas.scal(
                if (uplo == .upper)
                    j + 1
                else
                    n - j,
                beta,
                c + if (uplo == .upper)
                    j * ldc
                else
                    j + j * ldc,
                1,
            ) catch unreachable;
        }
    }

    if (transb == .no_trans or transb == .conj_no_trans) {
        if (transa == .no_trans or transa == .conj_no_trans) {
            for (0..n) |j| {
                const start = if (uplo == .upper) 0 else j;
                const end = if (uplo == .upper) j + 1 else n;
                for (0..k) |l| {
                    // temp = alpha * b[l + j * ldb]
                    const temp = numeric.mul(
                        alpha,
                        if (comptime noconjb)
                            b[l + j * ldb]
                        else
                            numeric.conj(b[l + j * ldb]),
                    );

                    for (start..end) |i| {
                        // c[i + j * ldc] += temp * a[i + l * lda]
                        numeric.fmaInto(
                            &c[i + j * ldc],
                            temp,
                            if (comptime noconja)
                                a[i + l * lda]
                            else
                                numeric.conj(a[i + l * lda]),
                            c[i + j * ldc],
                        );
                    }
                }
            }
        } else {
            for (0..n) |j| {
                const start = if (uplo == .upper) 0 else j;
                const end = if (uplo == .upper) j + 1 else n;
                for (start..end) |i| {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, B)));

                    for (0..k) |l| {
                        // temp += a[l + i * lda] * b[l + j * ldb]
                        numeric.fmaInto(
                            &temp,
                            if (comptime noconja)
                                a[l + i * lda]
                            else
                                numeric.conj(a[l + i * lda]),
                            if (comptime noconjb)
                                b[l + j * ldb]
                            else
                                numeric.conj(b[l + j * ldb]),
                            temp,
                        );
                    }

                    // c[i + j * ldc] += alpha * temp
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        alpha,
                        temp,
                        c[i + j * ldc],
                    );
                }
            }
        }
    } else {
        if (transa == .no_trans or transa == .conj_no_trans) {
            for (0..n) |j| {
                const start = if (uplo == .upper) 0 else j;
                const end = if (uplo == .upper) j + 1 else n;
                for (0..k) |l| {
                    // temp = alpha * b[j + l * ldb]
                    const temp = numeric.mul(
                        alpha,
                        if (comptime noconjb)
                            b[j + l * ldb]
                        else
                            numeric.conj(b[j + l * ldb]),
                    );

                    for (start..end) |i| {
                        // c[i + j * ldc] += temp * a[i + l * lda]
                        numeric.fmaInto(
                            &c[i + j * ldc],
                            if (comptime noconja)
                                a[i + l * lda]
                            else
                                numeric.conj(a[i + l * lda]),
                            temp,
                            c[i + j * ldc],
                        );
                    }
                }
            }
        } else {
            for (0..n) |j| {
                const start = if (uplo == .upper) 0 else j;
                const end = if (uplo == .upper) j + 1 else n;
                for (start..end) |i| {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, B)));

                    for (0..k) |l| {
                        // temp += a[l + i * lda] * b[j + l * ldb]
                        numeric.fmaInto(
                            &temp,
                            if (comptime noconja)
                                a[l + i * lda]
                            else
                                numeric.conj(a[l + i * lda]),
                            if (comptime noconjb)
                                b[j + l * ldb]
                            else
                                numeric.conj(b[j + l * ldb]),
                            temp,
                        );
                    }

                    // c[i + j * ldc] += alpha * temp
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        alpha,
                        temp,
                        c[i + j * ldc],
                    );
                }
            }
        }
    }

    return;
}
