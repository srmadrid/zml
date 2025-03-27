const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");

const Scalar = core.types.Scalar;

pub inline fn nrm2(comptime T: type, n: isize, x: [*]const T, incx: isize) Scalar(T) {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    if (n <= 0) return 0;

    const huge = std.math.floatMax(Scalar(T));
    var tsml: Scalar(T) = undefined;
    var tbig: Scalar(T) = undefined;
    var ssml: Scalar(T) = undefined;
    var sbig: Scalar(T) = undefined;
    if (Scalar(T) == f128) {
        tsml = core.math.pow128(2, @ceil((std.math.floatExponentMin(Scalar(T)) - 1) * 0.5));
        tbig = core.math.pow128(2, @floor((std.math.floatExponentMax(Scalar(T)) - @bitSizeOf(Scalar(T)) + 1) * 0.5));
        ssml = core.math.pow128(2, -@floor((std.math.floatExponentMin(Scalar(T)) - @bitSizeOf(Scalar(T))) * 0.5));
        sbig = core.math.pow128(2, -@ceil((std.math.floatExponentMax(Scalar(T)) + @bitSizeOf(Scalar(T)) - 1) * 0.5));
    } else {
        tsml = std.math.pow(Scalar(T), 2, @ceil((std.math.floatExponentMin(Scalar(T)) - 1) * 0.5));
        tbig = std.math.pow(Scalar(T), 2, @floor((std.math.floatExponentMax(Scalar(T)) - @bitSizeOf(Scalar(T)) + 1) * 0.5));
        ssml = std.math.pow(Scalar(T), 2, -@floor((std.math.floatExponentMin(Scalar(T)) - @bitSizeOf(Scalar(T))) * 0.5));
        sbig = std.math.pow(Scalar(T), 2, -@ceil((std.math.floatExponentMax(Scalar(T)) + @bitSizeOf(Scalar(T)) - 1) * 0.5));
    }

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
                    if (Scalar(T) == f128) {
                        abig += core.math.pow128(ax * sbig, 2);
                    } else {
                        abig += std.math.pow(Scalar(T), ax * sbig, 2);
                    }

                    notbig = false;
                } else if (ax < tsml) {
                    if (Scalar(T) == f128) {
                        if (notbig) asml += core.math.pow128(ax * ssml, 2);
                    } else {
                        if (notbig) asml += std.math.pow(Scalar(T), ax * ssml, 2);
                    }
                } else {
                    if (Scalar(T) == f128) {
                        amed += core.math.pow128(ax, 2);
                    } else {
                        amed += std.math.pow(Scalar(T), ax, 2);
                    }
                }
                ix += incx;
            }

            if (abig > 0) {
                if (amed > 0 or amed > huge or amed != amed) {
                    if (Scalar(T) == f128) {
                        abig += core.math.pow128(amed * sbig, 2);
                    } else {
                        abig += std.math.pow(Scalar(T), amed * sbig, 2);
                    }
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
                    if (Scalar(T) == f128) {
                        sumsq = core.math.pow128(ymax, 2) * (1 + core.math.pow128(ymin / ymax, 2));
                    } else {
                        sumsq = std.math.pow(Scalar(T), ymax, 2) * (1 + std.math.pow(Scalar(T), ymin / ymax, 2));
                    }
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
                    if (Scalar(T) == f128) {
                        abig += core.math.pow128(ax * sbig, 2);
                    } else {
                        abig += std.math.pow(Scalar(T), ax * sbig, 2);
                    }
                    notbig = false;
                } else if (ax < tsml) {
                    if (Scalar(T) == f128) {
                        if (notbig) asml += core.math.pow128(ax * ssml, 2);
                    } else {
                        if (notbig) asml += std.math.pow(Scalar(T), ax * ssml, 2);
                    }
                } else {
                    if (Scalar(T) == f128) {
                        amed += core.math.pow128(ax, 2);
                    } else {
                        amed += std.math.pow(Scalar(T), ax, 2);
                    }
                }
                ax = @abs(x[@intCast(ix)].im);
                if (ax > tbig) {
                    if (Scalar(T) == f128) {
                        abig += core.math.pow128(ax * sbig, 2);
                    } else {
                        abig += std.math.pow(Scalar(T), ax * sbig, 2);
                    }
                    notbig = false;
                } else if (ax < tsml) {
                    if (Scalar(T) == f128) {
                        if (notbig) asml += core.math.pow128(ax * ssml, 2);
                    } else {
                        if (notbig) asml += std.math.pow(Scalar(T), ax * ssml, 2);
                    }
                } else {
                    if (Scalar(T) == f128) {
                        amed += core.math.pow128(ax, 2);
                    } else {
                        amed += std.math.pow(Scalar(T), ax, 2);
                    }
                }
                ix += incx;
            }

            if (abig > 0) {
                if (amed > 0 or amed > huge or amed != amed) {
                    if (Scalar(T) == f128) {
                        abig += core.math.pow128(amed * sbig, 2);
                    } else {
                        abig += std.math.pow(Scalar(T), amed * sbig, 2);
                    }
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
                    if (Scalar(T) == f128) {
                        sumsq = core.math.pow128(ymax, 2) * (1 + core.math.pow128(ymin / ymax, 2));
                    } else {
                        sumsq = std.math.pow(Scalar(T), ymax, 2) * (1 + std.math.pow(Scalar(T), ymin / ymax, 2));
                    }
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
        .unsupported => unreachable,
    }

    return scl * @sqrt(sumsq);
}

test nrm2 {
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
