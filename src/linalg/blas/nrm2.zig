const std = @import("std");
const core = @import("../../core/core.zig");
const blas = @import("../blas.zig");

const scalar = core.supported.scalar;

pub inline fn nrm2(comptime T: type, n: isize, x: [*]const T, incx: isize) scalar(T) {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return 0;

    const huge = std.math.floatMax(scalar(T));
    const tsml = std.math.pow(scalar(T), 2, @ceil((std.math.floatExponentMin(f64) - 1) * 0.5));
    const tbig = std.math.pow(scalar(T), 2, @floor((std.math.floatExponentMax(f64) - @bitSizeOf(f64) + 1) * 0.5));
    const ssml = std.math.pow(scalar(T), 2, -@floor((std.math.floatExponentMin(f64) - @bitSizeOf(f64)) * 0.5));
    const sbig = std.math.pow(scalar(T), 2, -@ceil((std.math.floatExponentMax(f64) + @bitSizeOf(f64) - 1) * 0.5));

    var scl: scalar(T) = 1;
    var sumsq: scalar(T) = 0;

    var abig: scalar(T) = 0;
    var amed: scalar(T) = 0;
    var asml: scalar(T) = 0;
    var ix: isize = if (incx < 0) (1 - (n - 1) * incx) else 0;

    var notbig = true;

    switch (supported) {
        .BuiltinBool => @compileError("blas.nrm2 does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            var ax: scalar(T) = undefined;
            for (0..@intCast(n)) |_| {
                ax = @abs(x[@intCast(ix)]);
                if (ax > tbig) {
                    abig += std.math.pow(scalar(T), ax * sbig, 2);
                    notbig = false;
                } else if (ax < tsml) {
                    if (notbig) asml += std.math.pow(scalar(T), ax * ssml, 2);
                } else {
                    amed += std.math.pow(scalar(T), ax, 2);
                }
                ix += incx;
            }

            if (abig > 0) {
                if (amed > 0 or amed > huge or amed != amed) {
                    abig += std.math.pow(scalar(T), amed * sbig, 2);
                }
                scl = 1 / sbig;
                sumsq = abig;
            } else if (asml > 0) {
                if (amed > 0 or amed > huge or amed != amed) {
                    const sqrt_amed = @sqrt(amed);
                    const sqrt_asml = @sqrt(asml) / ssml;
                    const ymin = if (sqrt_asml > sqrt_amed) sqrt_amed else sqrt_asml;
                    const ymax = if (sqrt_asml > sqrt_amed) sqrt_asml else sqrt_amed;
                    scl = 1;
                    sumsq = std.math.pow(scalar(T), ymax, 2) * (1 + std.math.pow(scalar(T), ymin / ymax, 2));
                } else {
                    scl = 1 / ssml;
                    sumsq = asml;
                }
            } else {
                scl = 1;
                sumsq = amed;
            }
        },
        .Complex => {
            var ax: scalar(T) = undefined;
            for (0..@intCast(n)) |_| {
                ax = @abs(x[@intCast(ix)].re);
                if (ax > tbig) {
                    abig += std.math.pow(scalar(T), ax * sbig, 2);
                    notbig = false;
                } else if (ax < tsml) {
                    if (notbig) asml += std.math.pow(scalar(T), ax * ssml, 2);
                } else {
                    amed += std.math.pow(scalar(T), ax, 2);
                }
                ax = @abs(x[@intCast(ix)].im);
                if (ax > tbig) {
                    abig += std.math.pow(scalar(T), ax * sbig, 2);
                    notbig = false;
                } else if (ax < tsml) {
                    if (notbig) asml += std.math.pow(scalar(T), ax * ssml, 2);
                } else {
                    amed += std.math.pow(scalar(T), ax, 2);
                }
                ix += incx;
            }

            if (abig > 0) {
                if (amed > 0 or amed > huge or amed != amed) {
                    abig += std.math.pow(scalar(T), amed * sbig, 2);
                }
                scl = 1 / sbig;
                sumsq = abig;
            } else if (asml > 0) {
                if (amed > 0 or amed > huge or amed != amed) {
                    const sqrt_amed = @sqrt(amed);
                    const sqrt_asml = @sqrt(asml) / ssml;
                    const ymin = if (sqrt_asml > sqrt_amed) sqrt_amed else sqrt_asml;
                    const ymax = if (sqrt_asml > sqrt_amed) sqrt_asml else sqrt_amed;
                    scl = 1;
                    sumsq = std.math.pow(scalar(T), ymax, 2) * (1 + std.math.pow(scalar(T), ymin / ymax, 2));
                } else {
                    scl = 1 / ssml;
                    sumsq = asml;
                }
            } else {
                scl = 1;
                sumsq = amed;
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("blas.nrm2 only supports simple types."),
        .Unsupported => unreachable,
    }

    return scl * @sqrt(sumsq);
}

test "nrm2" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);

    for (0..n) |i| {
        x1[i] = @floatFromInt(i + 1);
    }

    const result1 = blas.nrm2(f64, n, x1.ptr, 1);
    try std.testing.expectApproxEqRel(18271.111077326415, result1, 0.0000000001);
    const result2 = blas.nrm2(f64, n / 2, x1.ptr, 2);
    try std.testing.expectApproxEqRel(12909.9380323842, result2, 0.0000000001);

    var x2 = try a.alloc(Complex(f64), n);
    defer a.free(x2);

    for (0..n) |i| {
        x2[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
    }

    const result3 = blas.nrm2(Complex(f64), n, x2.ptr, 1);
    try std.testing.expectApproxEqRel(25839.253085180306, result3, 0.0000000001);
    const result4 = blas.nrm2(Complex(f64), n / 2, x2.ptr, 2);
    try std.testing.expectApproxEqRel(18257.409454793964, result4, 0.0000000001);
}
