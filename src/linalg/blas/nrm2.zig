const std = @import("std");
const types = @import("../../types.zig");
const float = @import("../../float.zig");
const blas = @import("../blas.zig");

const Scalar = types.Scalar;

pub inline fn nrm2(comptime T: type, n: isize, x: [*]const T, incx: isize) Scalar(T) {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return 0;

    const huge = std.math.floatMax(Scalar(T));
    const tsml: Scalar(T) = float.pow(2, float.ceil((std.math.floatExponentMin(Scalar(T)) - 1) * @as(Scalar(T), 0.5)));
    const tbig: Scalar(T) = float.pow(2, float.floor((std.math.floatExponentMax(Scalar(T)) - @bitSizeOf(Scalar(T)) + 1) * @as(Scalar(T), 0.5)));
    const ssml: Scalar(T) = float.pow(2, -float.floor((std.math.floatExponentMin(Scalar(T)) - @bitSizeOf(Scalar(T))) * @as(Scalar(T), 0.5)));
    const sbig: Scalar(T) = float.pow(2, -float.ceil((std.math.floatExponentMax(Scalar(T)) + @bitSizeOf(Scalar(T)) - 1) * @as(Scalar(T), 0.5)));

    var scl: Scalar(T) = 1;
    var sumsq: Scalar(T) = 0;

    var abig: Scalar(T) = 0;
    var amed: Scalar(T) = 0;
    var asml: Scalar(T) = 0;
    var ix: isize = if (incx < 0) (1 - (n - 1) * incx) else 0;

    var notbig = true;

    switch (numericType) {
        .bool => @compileError("blas.nrm2 does not support bool."),
        .int, .float => {
            var ax: Scalar(T) = undefined;
            for (0..@intCast(n)) |_| {
                ax = @abs(x[@intCast(ix)]);
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
        },
        .cfloat => {
            var ax: Scalar(T) = undefined;
            for (0..@intCast(n)) |_| {
                ax = @abs(x[@intCast(ix)].re);
                if (ax > tbig) {
                    abig += float.pow(ax * sbig, 2);
                    notbig = false;
                } else if (ax < tsml) {
                    if (notbig) asml += float.pow(ax * ssml, 2);
                } else {
                    amed += float.pow(ax, 2);
                }
                ax = @abs(x[@intCast(ix)].im);
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
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.nrm2 only supports simple types."),
    }

    return scl * float.sqrt(sumsq);
}
