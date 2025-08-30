const std = @import("std");

const types = @import("../../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
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
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));
    const R: type = EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y));

    if (x.rows != y.rows or x.cols != y.cols)
        return matrix.Error.DimensionMismatch;

    var result: R = try .init(allocator, x.rows, x.cols);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(x)) == .col_major) { // orderOf(result) == orderOf(x)
        if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // c c c
            var j: u32 = 0;
            while (j < x.cols) : (j += 1) {
                var i: u32 = 0;
                while (i < if (j < x.upper) 0 else j - x.upper) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                    }
                }

                while (i <= int.min(x.rows - 1, j + x.lower)) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(x.data[(x.upper + i - j) + j * x.ld], y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(x.data[(x.upper + i - j) + j * x.ld], y.data[i + j * y.ld], ctx);
                    }
                }

                while (i < x.rows) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                    }
                }
            }
        } else { // c c r
            var j: u32 = 0;
            while (j < x.cols) : (j += 1) {
                var i: u32 = 0;
                while (i < if (j < x.upper) 0 else j - x.upper) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                    }
                }

                while (i <= int.min(x.rows - 1, j + x.lower)) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(x.data[(x.upper + i - j) + j * x.ld], y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(x.data[(x.upper + i - j) + j * x.ld], y.data[i * y.ld + j], ctx);
                    }
                }

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
        if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // r r c
            var i: u32 = 0;
            while (i < x.rows) : (i += 1) {
                var j: u32 = 0;
                while (j < if (i < x.lower) 0 else i - x.lower) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                    }
                }

                while (j <= int.min(x.cols - 1, i + x.upper)) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(x.data[i * x.ld + (x.lower + j - i)], y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + (x.lower + j - i)], y.data[i + j * y.ld], ctx);
                    }
                }

                while (j < x.cols) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                    }
                }
            }
        } else { // r r r
            var i: u32 = 0;
            while (i < x.rows) : (i += 1) {
                var j: u32 = 0;
                while (j < if (i < x.lower) 0 else i - x.lower) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                    }
                }

                while (j <= int.min(x.cols - 1, i + x.upper)) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(x.data[i * x.ld + (x.lower + j - i)], y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + (x.lower + j - i)], y.data[i * y.ld + j], ctx);
                    }
                }

                while (j < x.cols) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                    }
                }
            }
        }
    }

    return result;
}
