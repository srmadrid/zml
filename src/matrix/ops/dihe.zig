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

    var result: R = if (comptime types.isHermitianMatrix(R))
        try .init( // Hermitian
            allocator,
            x.rows,
        )
    else
        try .init( // General
            allocator,
            x.rows,
            x.cols,
        );
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // orderOf(result) == orderOf(y)
        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu d cu
            var j: u32 = 0;
            while (j < x.cols) : (j += 1) {
                var i: u32 = 0;
                while (i < j) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                    }

                    if (comptime types.isComplex(X)) { // Result is a general matrix
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[j + i * result.ld] = op(
                                constants.zero(X, ctx) catch unreachable,
                                ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                            );
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[j + i * result.ld] = try op(
                                constants.zero(X, ctx) catch unreachable,
                                ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                                ctx,
                            );
                        }
                    }
                }

                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[j + j * result.ld] = op(x.data[j], y.data[j + j * y.ld]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[j + j * result.ld] = try op(x.data[j], y.data[j + j * y.ld], ctx);
                }
            }
        } else { // cl d cl
            var j: u32 = 0;
            while (j < x.cols) : (j += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[j + j * result.ld] = op(x.data[j], y.data[j + j * y.ld]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[j + j * result.ld] = try op(x.data[j], y.data[j + j * y.ld], ctx);
                }

                var i: u32 = j + 1;
                while (i < x.rows) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                    }

                    if (comptime types.isComplex(X)) { // Result is a general matrix
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[j + i * result.ld] = op(
                                constants.zero(X, ctx) catch unreachable,
                                ops.conjugate(y.data[i + j * y.ld], ctx) catch unreachable,
                            );
                        } else if (comptime opinfo.@"fn".params.len == 3) {
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
        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru d ru
            var i: u32 = 0;
            while (i < x.rows) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i * result.ld + i] = op(x.data[i], y.data[i * y.ld + i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i * result.ld + i] = try op(x.data[i], y.data[i * y.ld + i], ctx);
                }

                var j: u32 = i + 1;
                while (j < x.cols) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                    }

                    if (comptime types.isComplex(X)) { // Result is a general matrix
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[j * result.ld + i] = op(
                                constants.zero(X, ctx) catch unreachable,
                                ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                            );
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[j * result.ld + i] = try op(
                                constants.zero(X, ctx) catch unreachable,
                                ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                                ctx,
                            );
                        }
                    }
                }
            }
        } else { // rl d rl
            var i: u32 = 0;
            while (i < x.rows) : (i += 1) {
                var j: u32 = 0;
                while (j < i) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                    }

                    if (comptime types.isComplex(X)) { // Result is a general matrix
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[j * result.ld + i] = op(
                                constants.zero(X, ctx) catch unreachable,
                                ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                            );
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[j * result.ld + i] = try op(
                                constants.zero(X, ctx) catch unreachable,
                                ops.conjugate(y.data[i * y.ld + j], ctx) catch unreachable,
                                ctx,
                            );
                        }
                    }
                }

                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i * result.ld + i] = op(x.data[i], y.data[i * y.ld + i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i * result.ld + i] = try op(x.data[i], y.data[i * y.ld + i], ctx);
                }
            }
        }
    }

    return result;
}
