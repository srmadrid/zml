const std = @import("std");

const types = @import("../types.zig");
const Coerce = types.Coerce;
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

// (x * y) + z computed with less rounding error than naively doing the
// multiplication and addition separately.
pub fn fma(x: anytype, y: anytype, z: anytype) Coerce(@TypeOf(x), Coerce(@TypeOf(y), @TypeOf(z))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Z: type = @TypeOf(z);
    const C: type = Coerce(X, Coerce(Y, Z));
    const S: type = types.Scalar(C);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float and types.numericType(X) != .cfloat) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float and types.numericType(Y) != .cfloat) or
        (types.numericType(Z) != .bool and types.numericType(Z) != .int and types.numericType(Z) != .float and types.numericType(Z) != .cfloat) or
        (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat and types.numericType(Z) != .cfloat))
        @compileError("cfloat.fma requires at least one of x, y or z to be a cfloat, the others must be bool, int, float or cfloat, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (comptime !types.isComplex(X)) {
        if (comptime !types.isComplex(Y)) {
            // z here is always complex
            return .{
                .re = float.fma(x, y, z.re),
                .im = types.scast(S, z.im),
            };
        } else {
            if (comptime !types.isComplex(Z)) {
                return .{
                    .re = float.fma(x, y.re, z),
                    .im = types.scast(S, x) * types.scast(S, y.im),
                };
            } else {
                return .{
                    .re = float.fma(x, y.re, z.re),
                    .im = float.fma(x, y.im, z.im),
                };
            }
        }
    } else {
        if (comptime !types.isComplex(Y)) {
            if (comptime !types.isComplex(Z)) {
                return .{
                    .re = float.fma(x.re, y, z),
                    .im = types.scast(S, x.im) * types.scast(S, y),
                };
            } else {
                return .{
                    .re = float.fma(x.re, y, z.re),
                    .im = float.fma(x.im, y, z.im),
                };
            }
        } else {
            if (comptime !types.isComplex(Z)) {
                return .{
                    .re = float.fma(x.re, y.re, float.fma(-x.im, y.im, z)),
                    .im = float.fma(x.re, y.im, x.im * y.re),
                };
            } else {
                return .{
                    .re = float.fma(x.re, y.re, float.fma(-x.im, y.im, z.re)),
                    .im = float.fma(x.re, y.im, float.fma(x.im, y.re, z.im)),
                };
            }
        }
    }
}
