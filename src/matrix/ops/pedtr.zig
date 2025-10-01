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

    var result: R = try .init(allocator, x.size, x.size);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // orderOf(result) == orderOf(y)
        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cu
            var j: u32 = 0;
            while (j < x.size) : (j += 1) {
                var i: u32 = 0;
                while (i < j) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(x.at(i, j), y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(x.at(i, j), y.data[i + j * y.ld], ctx);
                    }
                }

                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[j + j * result.ld] = op(x.at(j, j), constants.one(Y, ctx) catch unreachable);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[j + j * result.ld] = try op(x.at(j, j), constants.one(Y, ctx) catch unreachable, ctx);
                    }
                } else {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[j + j * result.ld] = op(x.at(j, j), y.data[j + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[j + j * result.ld] = try op(x.at(j, j), y.data[j + j * y.ld], ctx);
                    }
                }

                i = j + 1;
                while (i < x.size) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(x.at(i, j), constants.zero(Y, ctx) catch unreachable);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(x.at(i, j), constants.zero(Y, ctx) catch unreachable, ctx);
                    }
                }
            }
        } else { // c cl
            var j: u32 = 0;
            while (j < x.size) : (j += 1) {
                var i: u32 = 0;
                while (i < j) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(x.at(i, j), constants.zero(Y, ctx) catch unreachable);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(x.at(i, j), constants.zero(Y, ctx) catch unreachable, ctx);
                    }
                }

                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[j + j * result.ld] = op(x.at(j, j), constants.one(Y, ctx) catch unreachable);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[j + j * result.ld] = try op(x.at(j, j), constants.one(Y, ctx) catch unreachable, ctx);
                    }
                } else {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[j + j * result.ld] = op(x.at(j, j), y.data[j + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[j + j * result.ld] = try op(x.at(j, j), y.data[j + j * y.ld], ctx);
                    }
                }

                i = j + 1;
                while (i < x.size) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(x.at(i, j), y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(x.at(i, j), y.data[i + j * y.ld], ctx);
                    }
                }
            }
        }
    } else {
        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r ru
            var i: u32 = 0;
            while (i < x.size) : (i += 1) {
                var j: u32 = 0;
                while (j < i) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(x.at(i, j), constants.zero(Y, ctx) catch unreachable);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(x.at(i, j), constants.zero(Y, ctx) catch unreachable, ctx);
                    }
                }

                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + i] = op(x.at(i, i), constants.one(Y, ctx) catch unreachable);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + i] = try op(x.at(i, i), constants.one(Y, ctx) catch unreachable, ctx);
                    }
                } else {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + i] = op(x.at(i, i), y.data[i * y.ld + i]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + i] = try op(x.at(i, i), y.data[i * y.ld + i], ctx);
                    }
                }

                j = i + 1;
                while (j < x.size) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(x.at(i, j), y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(x.at(i, j), y.data[i * y.ld + j], ctx);
                    }
                }
            }
        } else { // r rl
            var i: u32 = 0;
            while (i < x.size) : (i += 1) {
                var j: u32 = 0;
                while (j < i) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(x.at(i, j), y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(x.at(i, j), y.data[i * y.ld + j], ctx);
                    }
                }

                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + i] = op(x.at(i, i), constants.one(Y, ctx) catch unreachable);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + i] = try op(x.at(i, i), constants.one(Y, ctx) catch unreachable, ctx);
                    }
                } else {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + i] = op(x.at(i, i), y.data[i * y.ld + i]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + i] = try op(x.at(i, i), y.data[i * y.ld + i], ctx);
                    }
                }

                j = i + 1;
                while (j < x.size) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(x.at(i, j), constants.zero(Y, ctx) catch unreachable);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(x.at(i, j), constants.zero(Y, ctx) catch unreachable, ctx);
                    }
                }
            }
        }
    }

    return result;
}
