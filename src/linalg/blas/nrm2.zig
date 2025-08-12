const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const Child = types.Child;
const EnsureFloat = types.EnsureFloat;
const float = @import("../../float.zig");
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn nrm2(
    n: i32,
    x: anytype,
    incx: i32,
    ctx: anytype,
) !EnsureFloat(Scalar(Child(@TypeOf(x)))) {
    const X: type = Child(@TypeOf(x));

    if (n <= 0) return blas.Error.InvalidArgument;

    if (comptime types.isArbitraryPrecision(X)) {
        // Orientative implementation for arbitrary precision types
        var sum: EnsureFloat(Scalar(X)) = try ops.init(EnsureFloat(Scalar(X)), ctx);
        errdefer ops.deinit(&sum, ctx);
        var temp: EnsureFloat(Scalar(X)) = try ops.init(EnsureFloat(Scalar(X)), ctx);
        defer ops.deinit(&temp, ctx);

        var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
        for (0..scast(u32, n)) |_| {
            if (comptime types.isComplex(X)) {
                try ops.abs2_(&temp, x[scast(u32, ix)], ctx);
            } else {
                try ops.pow_(&temp, x[scast(u32, ix)], 2, ctx);
            }

            try ops.add_(&sum, sum, temp, ctx);

            ix += incx;
        }

        try ops.sqrt_(&sum, sum, ctx);
        //return sum;

        @compileError("zml.linalg.blas.nrm2 not implemented for arbitrary precision types yet");
    } else {
        const huge = std.math.floatMax(EnsureFloat(Scalar(X)));
        const tsml: EnsureFloat(Scalar(X)) = float.pow(2, float.ceil((std.math.floatExponentMin(EnsureFloat(Scalar(X))) - 1) * @as(EnsureFloat(Scalar(X)), 0.5)));
        const tbig: EnsureFloat(Scalar(X)) = float.pow(2, float.floor((std.math.floatExponentMax(EnsureFloat(Scalar(X))) - @bitSizeOf(EnsureFloat(Scalar(X))) + 1) * @as(EnsureFloat(Scalar(X)), 0.5)));
        const ssml: EnsureFloat(Scalar(X)) = float.pow(2, -float.floor((std.math.floatExponentMin(EnsureFloat(Scalar(X))) - @bitSizeOf(EnsureFloat(Scalar(X)))) * @as(EnsureFloat(Scalar(X)), 0.5)));
        const sbig: EnsureFloat(Scalar(X)) = float.pow(2, -float.ceil((std.math.floatExponentMax(EnsureFloat(Scalar(X))) + @bitSizeOf(EnsureFloat(Scalar(X))) - 1) * @as(EnsureFloat(Scalar(X)), 0.5)));

        var scl: EnsureFloat(Scalar(X)) = 1;
        var sumsq: EnsureFloat(Scalar(X)) = 0;

        var abig: EnsureFloat(Scalar(X)) = 0;
        var amed: EnsureFloat(Scalar(X)) = 0;
        var asml: EnsureFloat(Scalar(X)) = 0;

        var notbig: bool = true;

        if (comptime types.numericType(X) == .cfloat) {
            var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
            for (0..scast(u32, n)) |_| {
                var ax: EnsureFloat(Scalar(X)) = float.abs(x[scast(u32, ix)].re);
                if (ax > tbig) {
                    abig += float.pow(ax * sbig, 2);
                    notbig = false;
                } else if (ax < tsml) {
                    if (notbig) asml += float.pow(ax * ssml, 2);
                } else {
                    amed += float.pow(ax, 2);
                }

                ax = float.abs(x[scast(u32, ix)].im);
                if (ax > tbig) {
                    abig += float.pow(ax * sbig, 2);
                    notbig = false;
                } else if (ax < tsml) {
                    if (notbig) asml += float.pow(ax * ssml, 2);
                } else {
                    amed += float.pow(ax, 2);
                }

                ix += incx;
            }
        } else {
            var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
            for (0..scast(u32, n)) |_| {
                const ax: EnsureFloat(Scalar(X)) = float.abs(x[scast(u32, ix)]);
                if (ax > tbig) {
                    abig += float.pow(ax * sbig, 2);
                    notbig = false;
                } else if (ax < tsml) {
                    if (notbig) asml += float.pow(ax * ssml, 2);
                } else {
                    amed += float.pow(ax, 2);
                }

                ix += incx;
            }
        }

        if (abig > 0) {
            if (amed > 0 or amed > huge or amed != amed) {
                abig += float.pow(amed * sbig, 2);
            }
            scl = 1 / sbig;
            sumsq = abig;
        } else if (asml > 0) {
            if (amed > 0 or amed > huge or amed != amed) {
                const sqrt_amed = float.sqrt(amed);
                const sqrt_asml = float.sqrt(asml) / ssml;
                const ymin = if (sqrt_asml > sqrt_amed) sqrt_amed else sqrt_asml;
                const ymax = if (sqrt_asml > sqrt_amed) sqrt_asml else sqrt_amed;
                scl = 1;
                sumsq = float.pow(ymax, 2) * (1 + float.pow(ymin / ymax, 2));
            } else {
                scl = 1 / ssml;
                sumsq = asml;
            }
        } else {
            scl = 1;
            sumsq = amed;
        }

        return scl * float.sqrt(sumsq);
    }
}
