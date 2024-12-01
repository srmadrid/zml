const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

pub inline fn scal(comptime T: type, n: isize, a: T, x: [*]T, incx: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) {
        return;
    }

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.scal does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (a == 1) return;

            if (incx == 1) {
                const m: usize = @as(usize, @intCast(n)) % 5;
                for (0..m) |i| {
                    x[i] *= a;
                }

                if (n < 5) return;

                var i: usize = m;
                while (i < n) : (i += 5) {
                    x[i] *= a;
                    x[i + 1] *= a;
                    x[i + 2] *= a;
                    x[i + 3] *= a;
                    x[i + 4] *= a;
                }
            } else {
                const nincx = n * incx;
                var i: isize = 0;
                while (i < nincx) : (i += incx) {
                    x[@intCast(i)] *= a;
                }
            }
        },
        .Complex => {
            if (a.re == 1 and a.im == 0) return;

            if (incx == 1) {
                const m: usize = @as(usize, @intCast(n)) % 5;
                for (0..m) |i| {
                    const temp = a.re * x[i].re - a.im * x[i].im;
                    x[i].im = a.re * x[i].im + a.im * x[i].re;
                    x[i].re = temp;
                }

                if (n < 5) return;

                var i: usize = m;
                while (i < n) : (i += 5) {
                    var temp = a.re * x[i].re - a.im * x[i].im;
                    x[i].im = a.re * x[i].im + a.im * x[i].re;
                    x[i].re = temp;
                    temp = a.re * x[i + 1].re - a.im * x[i + 1].im;
                    x[i + 1].im = a.re * x[i + 1].im + a.im * x[i + 1].re;
                    x[i + 1].re = temp;
                    temp = a.re * x[i + 2].re - a.im * x[i + 2].im;
                    x[i + 2].im = a.re * x[i + 2].im + a.im * x[i + 2].re;
                    x[i + 2].re = temp;
                    temp = a.re * x[i + 3].re - a.im * x[i + 3].im;
                    x[i + 3].im = a.re * x[i + 3].im + a.im * x[i + 3].re;
                    x[i + 3].re = temp;
                    temp = a.re * x[i + 4].re - a.im * x[i + 4].im;
                    x[i + 4].im = a.re * x[i + 4].im + a.im * x[i + 4].re;
                    x[i + 4].re = temp;
                }
            } else {
                const nincx = n * incx;
                var i: isize = 0;
                while (i < nincx) : (i += incx) {
                    const temp = a.re * x[@intCast(i)].re - a.im * x[@intCast(i)].im;
                    x[@intCast(i)].im = a.re * x[@intCast(i)].im + a.im * x[@intCast(i)].re;
                    x[@intCast(i)].re = temp;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.scal only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "scal" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 2, 3, 4 }, .{});
    defer A.deinit();

    A.setAll(1);

    try NDArray(f64).BLAS.scal(2, @constCast(&A.flatten()));

    for (0..A.size) |i| {
        try std.testing.expect(A.data[i] == 2);
    }

    const Complex = std.math.Complex;
    var B: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer B.deinit();

    B.setAll(Complex(f64).init(1, -1));

    try NDArray(Complex(f64)).BLAS.scal(Complex(f64).init(2, 1), @constCast(&B.flatten()));

    for (0..B.size) |i| {
        try std.testing.expect(std.math.approxEqAbs(f64, B.data[i].re, 3, 0.0001) and std.math.approxEqAbs(f64, B.data[i].im, -1, 0.0001));
    }
}
