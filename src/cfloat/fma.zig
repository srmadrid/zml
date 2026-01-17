const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const ops = @import("../ops.zig");

pub fn Fma(comptime X: type, comptime Y: type, comptime Z: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or !types.isNumeric(Z) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or !types.numericType(Z).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float and types.numericType(Z) != .float))
        @compileError("zml.cfloat.fma: at least one of x, y or z to be a cfloat, the others must be bool, int, float or cfloat, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tz: " ++ @typeName(Z) ++ "\n");

    return types.Coerce(X, types.Coerce(Y, Z));
}

// (x * y) + z computed with less rounding error than naively doing the
// multiplication and addition separately.
pub fn fma(x: anytype, y: anytype, z: anytype) Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Z: type = @TypeOf(z);
    const R: type = Fma(X, Y, Z);
    const S: type = types.Scalar(R);

    if (comptime !types.isComplex(X)) {
        if (comptime !types.isComplex(Y)) {
            // x, y are real, z is complex
            return .{
                .re = ops.fma(
                    x,
                    y,
                    z.re,
                    .{},
                ) catch unreachable,
                .im = types.scast(S, z.im),
            };
        } else {
            if (comptime !types.isComplex(Z)) {
                // x, z are real, y is complex
                return .{
                    .re = ops.fma(
                        x,
                        y.re,
                        z,
                        .{},
                    ) catch unreachable,
                    .im = ops.mul(
                        x,
                        y.im,
                        .{},
                    ) catch unreachable,
                };
            } else {
                // x is real, y, z are complex
                return .{
                    .re = ops.fma(
                        x,
                        y.re,
                        z.re,
                        .{},
                    ) catch unreachable,
                    .im = ops.fma(
                        x,
                        y.im,
                        z.im,
                        .{},
                    ) catch unreachable,
                };
            }
        }
    } else {
        if (comptime !types.isComplex(Y)) {
            if (comptime !types.isComplex(Z)) {
                // y, z are real, x is complex
                return .{
                    .re = ops.fma(
                        x.re,
                        y,
                        z,
                        .{},
                    ) catch unreachable,
                    .im = ops.mul(
                        x.im,
                        y,
                        .{},
                    ) catch unreachable,
                };
            } else {
                // y is real, x, z are complex
                return .{
                    .re = ops.fma(
                        x.re,
                        y,
                        z.re,
                        .{},
                    ) catch unreachable,
                    .im = ops.fma(
                        x.im,
                        y,
                        z.im,
                        .{},
                    ) catch unreachable,
                };
            }
        } else {
            if (comptime !types.isComplex(Z)) {
                // z is real, x, y are complex
                return .{
                    .re = ops.fma(
                        x.re,
                        y.re,
                        ops.fma(
                            -x.im,
                            y.im,
                            z,
                            .{},
                        ) catch unreachable,
                        .{},
                    ) catch unreachable,
                    .im = ops.fma(
                        x.re,
                        y.im,
                        ops.mul(
                            x.im,
                            y.re,
                            .{},
                        ) catch unreachable,
                        .{},
                    ) catch unreachable,
                };
            } else {
                // x, y, z are complex
                return .{
                    .re = ops.fma(
                        x.re,
                        y.re,
                        ops.fma(
                            ops.neg(x.im, .{}) catch unreachable,
                            y.im,
                            z.re,
                            .{},
                        ) catch unreachable,
                        .{},
                    ) catch unreachable,
                    .im = ops.fma(
                        x.re,
                        y.im,
                        ops.fma(
                            x.im,
                            y.re,
                            z.im,
                            .{},
                        ) catch unreachable,
                        .{},
                    ) catch unreachable,
                };
            }
        }
    }
}
