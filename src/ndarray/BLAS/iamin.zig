const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

pub inline fn iamin(comptime T: type, n: isize, x: [*]const T, incx: isize) usize {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0 or incx <= 0) return 0;

    if (n == 1) return 0;

    var imin: usize = 0;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.swap does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            var min = @abs(x[0]);

            if (incx == 1) {
                for (1..@intCast(n)) |i| {
                    if (@abs(x[i]) < @abs(min)) {
                        min = @abs(x[i]);
                        imin = i;
                    }
                }
            } else {
                var ix: isize = incx;

                for (1..@intCast(n)) |i| {
                    if (@abs(x[@intCast(ix)]) < @abs(min)) {
                        min = @abs(x[@intCast(ix)]);
                        imin = i;
                    }
                    ix += incx;
                }
            }
        },
        .Complex => {
            var min = @abs(x[0].re) + @abs(x[0].im);

            if (incx == 1) {
                for (1..@intCast(n)) |i| {
                    if (@abs(x[i].re) + @abs(x[i].im) < @abs(min)) {
                        min = @abs(x[i].re) + @abs(x[i].im);
                        imin = i;
                    }
                }
            } else {
                var ix: isize = incx;

                for (1..@intCast(n)) |i| {
                    if (@abs(x[@intCast(ix)].re) + @abs(x[@intCast(ix)].im) < @abs(min)) {
                        min = @abs(x[@intCast(ix)].re) + @abs(x[@intCast(ix)].im);
                        imin = i;
                    }
                    ix += incx;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.iamin only supports simple types."),
        .Unsupported => unreachable,
    }

    return imin;
}

test "iamin" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 2, 3, 4 }, .{});
    defer A.deinit();

    A.setAll(1);

    A.data[7] = 0;

    const result = try NDArray(f64).BLAS.iamin(A.flatten());

    try std.testing.expect(result == 7);

    const Complex = std.math.Complex;
    var C: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer C.deinit();

    C.setAll(Complex(f64).init(1, -1));

    C.data[7] = Complex(f64).init(1, 0);

    const result2 = try NDArray(Complex(f64)).BLAS.iamin(C.flatten());

    try std.testing.expect(result2 == 7);
}
