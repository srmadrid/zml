const std = @import("std");

const types = @import("../../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
const ops = @import("../../ops.zig");
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
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));
    const R: type = EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y));

    if (x.rows != y.size or x.cols != y.size)
        return matrix.Error.DimensionMismatch;

    var result: R = try .init(allocator, x.rows, x.cols);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(x)) == .col_major) { // orderOf(result) == orderOf(x)
        if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
            if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cu cu
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < j) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                result.data[j + i * result.ld] = op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                                );
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                result.data[j + i * result.ld] = try op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                                    ctx,
                                );
                            }
                        }

                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                            }
                        }
                    }
                } else { // c cu cl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < j) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(
                                    x.data[i + j * x.ld],
                                    ops.conjugate(y.data[j + i * y.ld], ctx) catch unreachable,
                                );
                                result.data[j + i * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[j + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(
                                    x.data[i + j * x.ld],
                                    ops.conjugate(y.data[j + i * y.ld], ctx) catch unreachable,
                                    ctx,
                                );
                                result.data[j + i * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[j + i * y.ld], ctx);
                            }
                        }

                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
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
            } else {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cl cu
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(
                                    x.data[i + j * x.ld],
                                    ops.conjugate(y.data[j + i * y.ld], ctx) catch unreachable,
                                );
                                result.data[j + i * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[j + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(
                                    x.data[i + j * x.ld],
                                    ops.conjugate(y.data[j + i * y.ld], ctx) catch unreachable,
                                    ctx,
                                );
                                result.data[j + i * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[j + i * y.ld], ctx);
                            }
                        }
                    }
                } else { // c cl cl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                result.data[j + i * result.ld] = op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                                );
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                result.data[j + i * result.ld] = try op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                                    ctx,
                                );
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cu ru
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < j) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                result.data[j + i * result.ld] = op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                                );
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                result.data[j + i * result.ld] = try op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                                    ctx,
                                );
                            }
                        }

                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                            }
                        }
                    }
                } else { // c cu rl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < j) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(
                                    x.data[i + j * x.ld],
                                    ops.conjugate(y.data[j * y.ld + i], ctx) catch unreachable,
                                );
                                result.data[j + i * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[j * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(
                                    x.data[i + j * x.ld],
                                    ops.conjugate(y.data[j * y.ld + i], ctx) catch unreachable,
                                    ctx,
                                );
                                result.data[j + i * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[j * y.ld + i], ctx);
                            }
                        }

                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
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
            } else {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cl ru
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(
                                    x.data[i + j * x.ld],
                                    ops.conjugate(y.data[j * y.ld + i], ctx) catch unreachable,
                                );
                                result.data[j + i * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[j * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(
                                    x.data[i + j * x.ld],
                                    ops.conjugate(y.data[j * y.ld + i], ctx) catch unreachable,
                                    ctx,
                                );
                                result.data[j + i * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[j * y.ld + i], ctx);
                            }
                        }
                    }
                } else { // c cl rl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                result.data[j + i * result.ld] = op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                                );
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                result.data[j + i * result.ld] = try op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                                    ctx,
                                );
                            }
                        }
                    }
                }
            }
        }
    } else {
        if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
            if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r ru cu
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                result.data[j * result.ld + i] = op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                                );
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                result.data[i + j * y.ld] = try op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                                    ctx,
                                );
                            }
                        }
                    }
                } else { // r ru cl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(
                                    x.data[i * x.ld + j],
                                    ops.conjugate(y.data[j + i * y.ld], ctx) catch unreachable,
                                );
                                result.data[j * result.ld + i] = op(constants.zero(X, ctx) catch unreachable, y.data[j + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(
                                    x.data[i * x.ld + j],
                                    ops.conjugate(y.data[j + i * y.ld], ctx) catch unreachable,
                                    ctx,
                                );
                                result.data[j * result.ld + i] = try op(constants.zero(X, ctx) catch unreachable, y.data[j + i * y.ld], ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r rl cu
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < i) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(
                                    x.data[i * x.ld + j],
                                    ops.conjugate(y.data[j + i * y.ld], ctx) catch unreachable,
                                );
                                result.data[j * result.ld + i] = op(constants.zero(X, ctx) catch unreachable, y.data[j + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(
                                    x.data[i * x.ld + j],
                                    ops.conjugate(y.data[j + i * y.ld], ctx) catch unreachable,
                                    ctx,
                                );
                                result.data[j * result.ld + i] = try op(constants.zero(X, ctx) catch unreachable, y.data[j + i * y.ld], ctx);
                            }
                        }

                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                            }
                        }
                    }
                } else { // r rl cl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < i) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                result.data[j * result.ld + i] = op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                                );
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                result.data[j * result.ld + i] = try op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                                    ctx,
                                );
                            }
                        }

                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
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
        } else {
            if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r ru ru
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                result.data[j * result.ld + i] = op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                                );
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                result.data[j * result.ld + i] = try op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                                    ctx,
                                );
                            }
                        }
                    }
                } else { // r ru rl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(
                                    x.data[i * x.ld + j],
                                    ops.conjugate(y.data[j * y.ld + i], ctx) catch unreachable,
                                );
                                result.data[j * result.ld + i] = op(constants.zero(X, ctx) catch unreachable, y.data[j * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(
                                    x.data[i * x.ld + j],
                                    ops.conjugate(y.data[j * y.ld + i], ctx) catch unreachable,
                                    ctx,
                                );
                                result.data[j * result.ld + i] = try op(constants.zero(X, ctx) catch unreachable, y.data[j * y.ld + i], ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r rl ru
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < i) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(
                                    x.data[i * x.ld + j],
                                    ops.conjugate(y.data[j * y.ld + i], ctx) catch unreachable,
                                );
                                result.data[j * result.ld + i] = op(constants.zero(X, ctx) catch unreachable, y.data[j * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(
                                    x.data[i * x.ld + j],
                                    ops.conjugate(y.data[j * y.ld + i], ctx) catch unreachable,
                                    ctx,
                                );
                                result.data[j * result.ld + i] = try op(constants.zero(X, ctx) catch unreachable, y.data[j * y.ld + i], ctx);
                            }
                        }

                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                            }
                        } else {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                            }
                        }
                    }
                } else { // r rl rl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < i) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                result.data[j * result.ld + i] = op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                                );
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                result.data[j * result.ld + i] = try op(
                                    constants.zero(X, ctx) catch unreachable,
                                    ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                                    ctx,
                                );
                            }
                        }

                        if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
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

    return result;
}
