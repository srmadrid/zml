const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

pub inline fn iamax(comptime T: type, n: isize, x: [*]const T, incx: isize) usize {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0 or incx <= 0) return 0;

    if (n == 1) return 0;

    var imax: usize = 0;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.iamax does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            var max = @abs(x[0]);

            if (incx == 1) {
                for (1..@intCast(n)) |i| {
                    if (@abs(x[i]) > @abs(max)) {
                        max = @abs(x[i]);
                        imax = i;
                    }
                }
            } else {
                var ix: isize = incx;

                for (1..@intCast(n)) |i| {
                    if (@abs(x[@intCast(ix)]) > @abs(max)) {
                        max = @abs(x[@intCast(ix)]);
                        imax = i;
                    }
                    ix += incx;
                }
            }
        },
        .Complex => {
            var max = @abs(x[0].re) + @abs(x[0].im);

            if (incx == 1) {
                for (1..@intCast(n)) |i| {
                    if (@abs(x[i].re) + @abs(x[i].im) > @abs(max)) {
                        max = @abs(x[i].re) + @abs(x[i].im);
                        imax = i;
                    }
                }
            } else {
                var ix: isize = incx;

                for (1..@intCast(n)) |i| {
                    if (@abs(x[@intCast(ix)].re) + @abs(x[@intCast(ix)].im) > @abs(max)) {
                        max = @abs(x[@intCast(ix)].re) + @abs(x[@intCast(ix)].im);
                        imax = i;
                    }
                    ix += incx;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.iamax only supports simple types."),
        .Unsupported => unreachable,
    }

    return imax;
}

test "iamax" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 2, 3, 4 }, .{});
    defer A.deinit();

    A.setAll(1);

    A.data[7] = 2;

    const result = try NDArray(f64).BLAS.iamax(A.flatten());

    try std.testing.expect(result == 7);

    const Complex = std.math.Complex;
    var C: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer C.deinit();

    C.setAll(Complex(f64).init(1, -1));

    C.data[7] = Complex(f64).init(2, 2);

    const result2 = try NDArray(Complex(f64)).BLAS.iamax(C.flatten());

    try std.testing.expect(result2 == 7);
}
