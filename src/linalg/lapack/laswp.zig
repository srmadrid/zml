const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;

pub inline fn laswp(
    order: Order,
    n: i32,
    a: anytype,
    lda: i32,
    k1: i32,
    k2: i32,
    ipiv: [*]const i32,
    incx: i32,
) !void {
    if (order == .col_major) {
        try k_laswp_c(n, a, lda, k1, k2, ipiv, incx);
    } else {
        try k_laswp_r(n, a, lda, k1, k2, ipiv, incx);
    }
}

fn k_laswp_c(
    n: i32,
    a: anytype,
    lda: i32,
    k1: i32,
    k2: i32,
    ipiv: [*]const i32,
    incx: i32,
) !void {
    const A: type = types.Child(@TypeOf(a));

    if (n <= 0 or incx == 0)
        return lapack.Error.InvalidArgument;

    // Interchange row i with row ipiv(k1 + (i - k1) * abs(incx)) for each of
    // rows k1 through k2.
    var ix0: i32 = undefined;
    var I1: i32 = undefined;
    var I2: i32 = undefined;
    var inc: i32 = undefined;
    if (incx > 0) {
        ix0 = k1 - 1;
        I1 = k1 - 1;
        I2 = k2 - 1;
        inc = 1;
    } else if (incx < 0) {
        ix0 = k1 - 1 + (k1 - k2) * incx;
        I1 = k2 - 1;
        I2 = k1 - 1;
        inc = -1;
    } else {
        return lapack.Error.InvalidArgument;
    }

    var temp: A = undefined;
    const n32: i32 = (n >> 5) << 5; // n32 = n - (n % 32)
    if (n32 != 0) {
        var j: i32 = 0;
        while (j < n32) {
            var ix: i32 = ix0;
            var i: i32 = I1;
            while (i != I2 + inc) {
                const ip: i32 = scast(i32, ipiv[scast(u32, ix)]) - 1; // Convert to 0-based index
                if (ip != i) {
                    var k: i32 = j;
                    while (k < j + 32) {
                        temp = a[scast(u32, i + k * lda)];
                        a[scast(u32, i + k * lda)] = a[scast(u32, ip + k * lda)];
                        a[scast(u32, ip + k * lda)] = temp;

                        k += 1;
                    }
                }

                i += inc;
                ix += incx;
            }

            j += 32;
        }
    }

    if (n32 != n) {
        var ix: i32 = ix0;
        var i: i32 = I1;
        while (i != I2 + inc) {
            const ip: i32 = scast(i32, ipiv[scast(u32, ix)]) - 1; // Convert to 0-based index
            if (ip != i) {
                var k: i32 = n32;
                while (k < n) {
                    temp = a[scast(u32, i + k * lda)];
                    a[scast(u32, i + k * lda)] = a[scast(u32, ip + k * lda)];
                    a[scast(u32, ip + k * lda)] = temp;

                    k += 1;
                }
            }

            i += inc;
            ix += incx;
        }
    }
}

fn k_laswp_r(
    n: i32,
    a: anytype,
    lda: i32,
    k1: i32,
    k2: i32,
    ipiv: [*]const i32,
    incx: i32,
) !void {
    const A: type = types.Child(@TypeOf(a));

    if (n <= 0 or lda < int.max(1, n) or incx == 0)
        return lapack.Error.InvalidArgument;

    // Interchange row i with row ipiv(k1 + (i - k1) * abs(incx)) for each of
    // rows k1 through k2.
    var ix0: i32 = undefined;
    var I1: i32 = undefined;
    var I2: i32 = undefined;
    var inc: i32 = undefined;
    if (incx > 0) {
        ix0 = k1 - 1;
        I1 = k1 - 1;
        I2 = k2 - 1;
        inc = 1;
    } else {
        ix0 = k1 - 1 + (k1 - k2) * incx;
        I1 = k2 - 1;
        I2 = k1 - 1;
        inc = -1;
    }

    var temp: A = undefined;
    const n32: i32 = (n >> 5) << 5; // n32 = n - (n % 32)
    if (n32 != 0) {
        var j: i32 = 0;
        while (j < n32) {
            var ix: i32 = ix0;
            var i: i32 = I1;
            while (i != I2 + inc) {
                const ip: i32 = scast(i32, ipiv[scast(u32, ix)]) - 1; // Convert to 0-based index
                if (ip != i) {
                    var k: i32 = j;
                    while (k < j + 32) {
                        temp = a[scast(u32, i * lda + k)];
                        a[scast(u32, i * lda + k)] = a[scast(u32, ip * lda + k)];
                        a[scast(u32, ip * lda + k)] = temp;

                        k += 1;
                    }
                }

                i += inc;
                ix += incx;
            }

            j += 32;
        }
    }

    if (n32 != n) {
        var ix: i32 = ix0;
        var i: i32 = I1;
        while (i != I2 + inc) {
            const ip: i32 = scast(i32, ipiv[scast(u32, ix)]) - 1; // Convert to 0-based index
            if (ip != i) {
                var k: i32 = n32;
                while (k < n) {
                    temp = a[scast(u32, i * lda + k)];
                    a[scast(u32, i * lda + k)] = a[scast(u32, ip * lda + k)];
                    a[scast(u32, ip * lda + k)] = temp;

                    k += 1;
                }
            }

            i += inc;
            ix += incx;
        }
    }
}
