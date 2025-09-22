const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Coerce = types.Coerce;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Mach = lapack.Mach;

pub fn lapy2(
    x: anytype,
    y: anytype,
    ctx: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    if (comptime !types.isArbitraryPrecision(C)) {
        if (types.numericType(C) == .int) {
            return types.scast(C, try ops.hypot(x, y, ctx));
        } else { // float
            if (std.math.isNan(x))
                return types.scast(C, x);

            if (std.math.isNan(y))
                return types.scast(C, y);

            const hugeval: C = lapack.lamch(C, .rmax);
            const xabs: C = try ops.abs(types.scast(C, x), ctx);
            const yabs: C = try ops.abs(types.scast(C, y), ctx);
            const w: C = try ops.max(xabs, yabs, ctx);
            const z: C = try ops.min(xabs, yabs, ctx);
            if (z == 0 or w > hugeval) {
                return w;
            } else {
                return w * try ops.sqrt(1 + try ops.pow(z / w, 2, ctx), ctx);
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.lapy2 not implemented for arbitrary precision types yet");
    }
}
