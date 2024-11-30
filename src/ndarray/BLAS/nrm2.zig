const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const Error = @import("../ndarray.zig").Error;
const core = @import("../../core/core.zig");

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
        .BuiltinBool => @compileError("BLAS.nrm2 does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            // Compute the sum of squares.
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

            // Combine the accumulators.
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
            // Compute the sum of squares.
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

            // Combine the accumulators.
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
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.nrm2 only supports simple types."),
        .Unsupported => unreachable,
    }

    return scl * @sqrt(sumsq);
}

test "nrm2" {
    const a = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 8, 18, 7 }, .{});
    defer A.deinit();

    A.setAll(1);

    const result1 = try NDArray(f64).BLAS.nrm2(A.flatten());
    try std.testing.expect(result1 == @sqrt(1008.0));

    const Complex = std.math.Complex;
    var B: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer B.deinit();

    B.setAll(Complex(f64).init(1, -1));

    const result2: f64 = try NDArray(Complex(f64)).BLAS.nrm2(B.flatten());
    try std.testing.expect(result2 == @sqrt(2016.0));
}
