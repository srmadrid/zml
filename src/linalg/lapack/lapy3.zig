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

pub fn lapy3(
    x: anytype,
    y: anytype,
    z: anytype,
    ctx: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Z: type = @TypeOf(z);
    const C: type = Coerce(X, types.Coerce(Y, Z));

    if (comptime !types.isArbitraryPrecision(C)) {
        if (types.numericType(C) == .int) {
            return try ops.sqrt(
                try ops.add(
                    try ops.add(
                        try ops.pow(types.scast(C, x), 2, ctx),
                        try ops.pow(types.scast(C, y), 2, ctx),
                        ctx,
                    ),
                    try ops.pow(types.scast(C, z), 2, ctx),
                    ctx,
                ),
                ctx,
            );
        } else { // float
            const hugeval: C = lapack.lamch(C, .rmax);
            const xabs: C = try ops.abs(types.scast(C, x), ctx);
            const yabs: C = try ops.abs(types.scast(C, y), ctx);
            const zabs: C = try ops.abs(types.scast(C, z), ctx);
            const w: C = try ops.max(xabs, try ops.max(yabs, zabs, ctx), ctx);
            if (w == 0 or w > hugeval) {
                return try ops.add(
                    try ops.add(
                        xabs,
                        yabs,
                        ctx,
                    ),
                    zabs,
                    ctx,
                );
            } else {
                return try ops.mul(
                    w,
                    try ops.sqrt(
                        try ops.add(
                            try ops.add(
                                try ops.mul(try ops.div(xabs, w, ctx), try ops.div(xabs, w, ctx), ctx),
                                try ops.mul(try ops.div(yabs, w, ctx), try ops.div(yabs, w, ctx), ctx),
                                ctx,
                            ),
                            try ops.mul(try ops.div(zabs, w, ctx), try ops.div(zabs, w, ctx), ctx),
                            ctx,
                        ),
                        ctx,
                    ),
                    ctx,
                );
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.lapy3 not implemented for arbitrary precision types yet");
    }
}
