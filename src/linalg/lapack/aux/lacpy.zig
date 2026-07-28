const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Uplo = types.Uplo;

const utils = @import("../utils.zig");

pub fn lacpy(
    order: Order,
    uplo: union(enum) { uplo: Uplo, full: void },
    m: i32,
    n: i32,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !void {
    switch (uplo) {
        .uplo => |_uplo| {
            if (order == .col_major) {
                if (_uplo == .upper) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = 0;
                        while (i <= int.min(j, m - 1)) : (i += 1) {
                            try ops.set(
                                &b[types.scast(u32, i + j * ldb)],
                                a[types.scast(u32, i + j * lda)],
                                ctx,
                            );
                        }
                    }
                } else {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = j;
                        while (i < m) : (i += 1) {
                            try ops.set(
                                &b[types.scast(u32, i + j * ldb)],
                                a[types.scast(u32, i + j * lda)],
                                ctx,
                            );
                        }
                    }
                }
            } else {
                if (_uplo == .upper) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        var j: i32 = i;
                        while (j < n) : (j += 1) {
                            try ops.set(
                                &b[types.scast(u32, i * ldb + j)],
                                a[types.scast(u32, i * lda + j)],
                                ctx,
                            );
                        }
                    }
                } else {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        var j: i32 = 0;
                        while (j <= int.min(i, n - 1)) : (j += 1) {
                            try ops.set(
                                &b[types.scast(u32, i * ldb + j)],
                                a[types.scast(u32, i * lda + j)],
                                ctx,
                            );
                        }
                    }
                }
            }
        },
        .full => {
            if (order == .col_major) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        try ops.set(
                            &b[types.scast(u32, i + j * ldb)],
                            a[types.scast(u32, i + j * lda)],
                            ctx,
                        );
                    }
                }
            } else {
                var i: i32 = 0;
                while (i < m) : (i += 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        try ops.set(
                            &b[types.scast(u32, i * ldb + j)],
                            a[types.scast(u32, i * lda + j)],
                            ctx,
                        );
                    }
                }
            }
        },
    }
}
