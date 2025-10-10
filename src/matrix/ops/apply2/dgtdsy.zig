const std = @import("std");

const types = @import("../../../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
const ops = @import("../../../ops.zig");
const constants = @import("../../../constants.zig");
const int = @import("../../../int.zig");

const matrix = @import("../../../matrix.zig");

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

    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, x.size, x.size);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // orderOf(result) == orderOf(y)
        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cu
            var j: u32 = 0;
            while (j < x.size) : (j += 1) {
                var i: u32 = 0;
                while (i < j) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(x.get(i, j) catch unreachable, y.data[i + j * y.ld]);
                        result.data[j + i * result.ld] = op(x.get(j, i) catch unreachable, y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(x.get(i, j) catch unreachable, y.data[i + j * y.ld], ctx);
                        result.data[j + i * result.ld] = try op(x.get(j, i) catch unreachable, y.data[i + j * y.ld], ctx);
                    }
                }

                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[j + j * result.ld] = op(x.get(j, j) catch unreachable, y.data[j + j * y.ld]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[j + j * result.ld] = try op(x.get(j, j) catch unreachable, y.data[j + j * y.ld], ctx);
                }
            }
        } else { // c cl
            var j: u32 = 0;
            while (j < x.size) : (j += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[j + j * result.ld] = op(x.get(j, j) catch unreachable, y.data[j + j * y.ld]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[j + j * result.ld] = try op(x.get(j, j) catch unreachable, y.data[j + j * y.ld], ctx);
                }

                var i: u32 = j + 1;
                while (i < x.size) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(x.get(i, j) catch unreachable, y.data[i + j * y.ld]);
                        result.data[j + i * result.ld] = op(x.get(j, i) catch unreachable, y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(x.get(i, j) catch unreachable, y.data[i + j * y.ld], ctx);
                        result.data[j + i * result.ld] = try op(x.get(j, i) catch unreachable, y.data[i + j * y.ld], ctx);
                    }
                }
            }
        }
    } else {
        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r ru
            var i: u32 = 0;
            while (i < x.size) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i * result.ld + i] = op(x.get(i, i) catch unreachable, y.data[i * y.ld + i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i * result.ld + i] = try op(x.get(i, i) catch unreachable, y.data[i * y.ld + i], ctx);
                }

                var j: u32 = i + 1;
                while (j < x.size) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(x.get(i, j) catch unreachable, y.data[i * y.ld + j]);
                        result.data[j * result.ld + i] = op(x.get(j, i) catch unreachable, y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(x.get(i, j) catch unreachable, y.data[i * y.ld + j], ctx);
                        result.data[j * result.ld + i] = try op(x.get(j, i) catch unreachable, y.data[i * y.ld + j], ctx);
                    }
                }
            }
        } else { // r rl
            var i: u32 = 0;
            while (i < x.size) : (i += 1) {
                var j: u32 = 0;
                while (j < i) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(x.get(i, j) catch unreachable, y.data[i * y.ld + j]);
                        result.data[j * result.ld + i] = op(x.get(j, i) catch unreachable, y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(x.get(i, j) catch unreachable, y.data[i * y.ld + j], ctx);
                        result.data[j * result.ld + i] = try op(x.get(j, i) catch unreachable, y.data[i * y.ld + j], ctx);
                    }
                }

                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i * result.ld + i] = op(x.get(i, i) catch unreachable, y.data[i * y.ld + i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i * result.ld + i] = try op(x.get(i, i) catch unreachable, y.data[i * y.ld + i], ctx);
                }
            }
        }
    }

    return result;
}
