const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

const scalar = core.supported.scalar;

pub inline fn rotm(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, param: [*]const T) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n == 0 or param[0] + 2 == 0) return;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.rotm does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            const flag = param[0];
            var kx: isize = 0;
            var ky: isize = 0;

            if (incx == incy and incx > 0) {
                const nsteps: usize = @intCast(n * incx);
                if (flag < 0) {
                    const h11 = param[1];
                    const h12 = param[3];
                    const h21 = param[2];
                    const h22 = param[4];

                    var i: usize = 0;
                    while (i < nsteps) : (i += @intCast(incx)) {
                        const idx = i * @as(usize, @intCast(incx));
                        const w = x[idx];
                        const z = y[idx];
                        x[idx] = w * h11 + z * h12;
                        y[idx] = w * h21 + z * h22;
                    }
                } else if (flag == 0.0) {
                    const h12 = param[3];
                    const h21 = param[2];

                    var i: usize = 0;
                    while (i < nsteps) : (i += @intCast(incx)) {
                        const idx = i * @as(usize, @intCast(incx));
                        const w = x[idx];
                        const z = y[idx];
                        x[idx] = w + z * h12;
                        y[idx] = w * h21 + z;
                    }
                } else {
                    const h11 = param[1];
                    const h22 = param[4];

                    var i: usize = 0;
                    while (i < nsteps) : (i += @intCast(incx)) {
                        const idx = i * @as(usize, @intCast(incx));
                        const w = x[idx];
                        const z = y[idx];
                        x[idx] = w * h11 + z;
                        y[idx] = -w + h22 * z;
                    }
                }
            } else {
                if (incx < 0) kx = (1 - n) * incx;
                if (incy < 0) ky = (1 - n) * incy;

                if (flag < 0) {
                    const h11 = param[1];
                    const h12 = param[3];
                    const h21 = param[2];
                    const h22 = param[4];

                    for (0..@intCast(n)) |_| {
                        const w = x[@intCast(kx)];
                        const z = y[@intCast(ky)];
                        x[@intCast(kx)] = w * h11 + z * h12;
                        y[@intCast(ky)] = w * h21 + z * h22;
                        kx += incx;
                        ky += incy;
                    }
                } else if (flag == 0) {
                    const h12 = param[3];
                    const h21 = param[2];

                    for (0..@intCast(n)) |_| {
                        const w = x[@intCast(kx)];
                        const z = y[@intCast(ky)];
                        x[@intCast(kx)] = w + z * h12;
                        y[@intCast(ky)] = w * h21 + z;
                        kx += incx;
                        ky += incy;
                    }
                } else {
                    const h11 = param[1];
                    const h22 = param[4];

                    for (0..@intCast(n)) |_| {
                        const w = x[@intCast(kx)];
                        const z = y[@intCast(ky)];
                        x[@intCast(kx)] = w * h11 + z;
                        y[@intCast(ky)] = -w + h22 * z;
                        kx += incx;
                        ky += incy;
                    }
                }
            }
        },
        .Complex => @compileError("BLAS.rotm does not support complex numbers."),
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.rotm only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "rotm" {
    const a = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 8, 18, 7 }, .{});
    defer A.deinit();

    A.setAll(1);

    var B: NDArray(f64) = try NDArray(f64).init(a, &.{ 8, 18, 7 }, .{});
    defer B.deinit();

    B.setAll(1);

    var param: NDArray(f64) = try NDArray(f64).init(a, &.{5}, .{});
    defer param.deinit();

    param.data[0] = -1;
    param.data[1] = 0.7071067811865475;
    param.data[2] = 0;
    param.data[3] = 0.7071067811865475;
    param.data[4] = 0;

    try NDArray(f64).BLAS.rotm(@constCast(&A.flatten()), @constCast(&B.flatten()), param.flatten());

    for (0..A.size) |i| {
        try std.testing.expectApproxEqAbs(1.414214, A.data[i], 0.0001);
        try std.testing.expectApproxEqAbs(0, B.data[i], 0.0001);
    }
}
