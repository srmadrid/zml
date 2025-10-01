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

    var result: R = try .init(allocator, x.size, x.size);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(x)) == .col_major) { // orderOf(result) == orderOf(x)
        if (comptime types.uploOf(@TypeOf(x)) == .upper) { // c cu
            var j: u32 = 0;
            while (j < x.size) : (j += 1) {
                var i: u32 = 0;
                while (i < j) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.get(i, j) catch unreachable);
                        result.data[j + i * result.ld] = op(
                            ops.conjugate(x.data[i + j * x.ld], ctx) catch unreachable,
                            y.get(j, i) catch unreachable,
                        );
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.get(i, j) catch unreachable, ctx);
                        result.data[j + i * result.ld] = try op(
                            ops.conjugate(x.data[i + j * x.ld], ctx) catch unreachable,
                            y.get(j, i) catch unreachable,
                            ctx,
                        );
                    }
                }

                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.get(j, j) catch unreachable);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.get(j, j) catch unreachable, ctx);
                }
            }
        } else { // c cl
            var j: u32 = 0;
            while (j < x.size) : (j += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.get(j, j) catch unreachable);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.get(j, j) catch unreachable, ctx);
                }

                var i: u32 = j + 1;
                while (i < x.size) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.get(i, j) catch unreachable);
                        result.data[j + i * result.ld] = op(
                            ops.conjugate(x.data[i + j * x.ld], ctx) catch unreachable,
                            y.get(j, i) catch unreachable,
                        );
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.get(i, j) catch unreachable, ctx);
                        result.data[j + i * result.ld] = try op(
                            ops.conjugate(x.data[i + j * x.ld], ctx) catch unreachable,
                            y.get(j, i) catch unreachable,
                            ctx,
                        );
                    }
                }
            }
        }
    } else {
        if (comptime types.uploOf(@TypeOf(x)) == .upper) { // r ru
            var i: u32 = 0;
            while (i < x.size) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.get(i, i) catch unreachable);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.get(i, i) catch unreachable, ctx);
                }

                var j: u32 = i + 1;
                while (j < x.size) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.get(i, j) catch unreachable);
                        result.data[j * result.ld + i] = op(
                            ops.conjugate(x.data[i * x.ld + j], ctx) catch unreachable,
                            y.get(j, i) catch unreachable,
                        );
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.get(i, j) catch unreachable, ctx);
                        result.data[j * result.ld + i] = try op(
                            ops.conjugate(x.data[i * x.ld + j], ctx) catch unreachable,
                            y.get(j, i) catch unreachable,
                            ctx,
                        );
                    }
                }
            }
        } else { // r rl
            var i: u32 = 0;
            while (i < x.size) : (i += 1) {
                var j: u32 = 0;
                while (j < i) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.get(i, j) catch unreachable);
                        result.data[j * result.ld + i] = op(
                            ops.conjugate(x.data[i * x.ld + j], ctx) catch unreachable,
                            y.get(j, i) catch unreachable,
                        );
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.get(i, j) catch unreachable, ctx);
                        result.data[j * result.ld + i] = try op(
                            ops.conjugate(x.data[i * x.ld + j], ctx) catch unreachable,
                            y.get(j, i) catch unreachable,
                            ctx,
                        );
                    }
                }

                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.get(i, i) catch unreachable);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.get(i, i) catch unreachable, ctx);
                }
            }
        }
    }

    return result;
}
