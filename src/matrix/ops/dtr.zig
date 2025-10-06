const std = @import("std");

const types = @import("../../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
const Uplo = types.Uplo;
const ops = @import("../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = ReturnType2(op, Numeric(X), Numeric(Y));

    if (comptime !types.isTriangularDenseMatrix(@TypeOf(x))) {
        var result: matrix.triangular.Dense(R, types.uploOf(Y), .non_unit, types.orderOf(Y)) = try .init(allocator, y.rows, y.cols);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cu
                    var j: u32 = 0;
                    while (j < y.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, y.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                            }
                        }

                        if (j < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                                }
                            }
                        }
                    }
                } else { // c cl
                    var j: u32 = 0;
                    while (j < y.cols) : (j += 1) {
                        if (j < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < y.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c ru
                    var j: u32 = 0;
                    while (j < y.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, y.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                            }
                        }

                        if (j < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                                }
                            }
                        }
                    }
                } else { // c rl
                    var j: u32 = 0;
                    while (j < y.cols) : (j += 1) {
                        if (j < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < y.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r cu
                    var i: u32 = 0;
                    while (i < y.rows) : (i += 1) {
                        if (i < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < y.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                            }
                        }
                    }
                } else { // r cl
                    var i: u32 = 0;
                    while (i < y.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, y.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                            }
                        }

                        if (i < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r ru
                    var i: u32 = 0;
                    while (i < y.rows) : (i += 1) {
                        if (i < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < y.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                            }
                        }
                    }
                } else { // r rl
                    var i: u32 = 0;
                    while (i < y.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, y.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                            }
                        }

                        if (i < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                                }
                            }
                        }
                    }
                }
            }
        }

        return result;
    } else if (comptime !types.isTriangularDenseMatrix(@TypeOf(y))) {
        var result: matrix.triangular.Dense(R, types.uploOf(X), .non_unit, types.orderOf(X)) = try .init(allocator, x.rows, x.cols);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // c cu
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, x.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                            }
                        }

                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                } else { // c cl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // c ru
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, x.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                            }
                        }

                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                                }
                            }
                        }
                    }
                } else { // c rl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // r cu
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                            }
                        }
                    }
                } else { // r cl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, x.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                            }
                        }

                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // r ru
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                            }
                        }
                    }
                } else { // r rl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, x.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                            }
                        }

                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
                                }
                            }
                        }
                    }
                }
            }
        }

        return result;
    }

    var result: if (types.uploOf(X) == types.uploOf(Y))
        matrix.triangular.Dense(R, types.uploOf(X), .non_unit, types.orderOf(Y))
    else
        matrix.general.Dense(R, types.orderOf(X)) =
        try .init(allocator, x.rows, x.cols);
    errdefer result.deinit(allocator);

    // Two cases:
    // - Result is triangular
    // - Result is general
    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.isTriangularDenseMatrix(@TypeOf(result))) { // uploOf(x) == uploOf(y)
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cu cu cu
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, x.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                            }
                        }

                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else { // cl cl cl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                                    }
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cu cu ru
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, x.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                            }
                        }

                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else { // cl cl rl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                                    }
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // ru ru cu
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                                    }
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                            }
                        }
                    }
                } else { // rl rl cl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, x.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                            }
                        }

                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // ru ru ru
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                                    }
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                            }
                        }
                    }
                } else { // rl rl rl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, x.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                            }
                        }

                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } else { // uploOf(x) != uploOf(y)
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cu cu
                        unreachable; // Handled above
                    } else { // c cu cl
                        var j: u32 = 0;
                        while (j < x.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, x.rows)) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }

                            if (j < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                                        }
                                    }
                                }
                            }

                            i = j + 1;
                            while (i < x.rows) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cl cu
                        var j: u32 = 0;
                        while (j < x.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, x.rows)) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                                }
                            }

                            if (j < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                                        }
                                    }
                                }
                            }

                            i = j + 1;
                            while (i < x.rows) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }
                        }
                    } else { // c cl cl
                        unreachable; // Handled above
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cu ru
                        unreachable; // Handled above
                    } else { // c cu rl
                        var j: u32 = 0;
                        while (j < x.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, x.rows)) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }

                            if (j < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                                        }
                                    }
                                }
                            }

                            i = j + 1;
                            while (i < x.rows) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cl ru
                        var j: u32 = 0;
                        while (j < x.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, x.rows)) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                                }
                            }

                            if (j < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                                        }
                                    }
                                }
                            }

                            i = j + 1;
                            while (i < x.rows) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }
                        }
                    } else { // c cl rl
                        unreachable; // Handled above
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r ru cu
                        unreachable; // Handled above
                    } else { // r ru cl
                        var i: u32 = 0;
                        while (i < x.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, x.cols)) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                                }
                            }

                            if (i < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                                        }
                                    }
                                }
                            }

                            j = i + 1;
                            while (j < x.cols) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r rl cu
                        var i: u32 = 0;
                        while (i < x.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, x.cols)) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }

                            if (i < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                                        }
                                    }
                                }
                            }

                            j = i + 1;
                            while (j < x.cols) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    } else { // r rl cl
                        unreachable; // Handled above
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r ru ru
                        unreachable; // Handled above
                    } else { // r ru rl
                        var i: u32 = 0;
                        while (i < x.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, x.cols)) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                                }
                            }

                            if (i < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                                        }
                                    }
                                }
                            }

                            j = i + 1;
                            while (j < x.cols) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r rl ru
                        var i: u32 = 0;
                        while (i < x.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, x.cols)) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }

                            if (i < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                                        }
                                    }
                                }
                            }

                            j = i + 1;
                            while (j < x.cols) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    } else { // r rl rl
                        unreachable; // Handled above
                    }
                }
            }
        }
    }

    return result;
}
