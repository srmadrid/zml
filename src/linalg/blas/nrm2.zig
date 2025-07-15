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
    n: isize,
    x: anytype,
    incx: isize,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !EnsureFloat(Scalar(Child(@TypeOf(x)))) {
    const X: type = Child(@TypeOf(x));

    if (n <= 0) return blas.Error.InvalidArgument;

    switch (comptime types.numericType(X)) {
        .int, .float, .cfloat => {
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
                var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
                for (0..scast(usize, n)) |_| {
                    var ax: EnsureFloat(Scalar(X)) = float.abs(x[scast(usize, ix)].re);
                    if (ax > tbig) {
                        abig += float.pow(ax * sbig, 2);
                        notbig = false;
                    } else if (ax < tsml) {
                        if (notbig) asml += float.pow(ax * ssml, 2);
                    } else {
                        amed += float.pow(ax, 2);
                    }

                    ax = float.abs(x[scast(usize, ix)].im);
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
                var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
                for (0..scast(usize, n)) |_| {
                    const ax: EnsureFloat(Scalar(X)) = float.abs(x[scast(usize, ix)]);
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
        },
        else => {
            var sum: EnsureFloat(Scalar(X)) = try ops.init(EnsureFloat(Scalar(X)), .{ .allocator = options.allocator });
            errdefer ops.deinit(&sum, .{ .allocator = options.allocator });
            var temp: EnsureFloat(Scalar(X)) = try ops.init(EnsureFloat(Scalar(X)), .{ .allocator = options.allocator });
            defer ops.deinit(&temp, .{ .allocator = options.allocator });

            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            for (0..scast(usize, n)) |_| {
                if (comptime types.isComplex(X)) {
                    try ops.abs2_(&temp, x[scast(usize, ix)], .{ .allocator = options.allocator });
                } else {
                    try ops.pow_(&temp, x[scast(usize, ix)], 2, .{ .allocator = options.allocator });
                }

                try ops.add_(&sum, sum, temp, .{ .allocator = options.allocator });

                ix += incx;
            }

            try ops.sqrt_(&sum, sum, .{ .allocator = options.allocator });
            return sum;
        },
    }
}
