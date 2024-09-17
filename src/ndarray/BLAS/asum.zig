const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const Error = @import("../ndarray.zig").Error;
const core = @import("../../core/core.zig");

const _add = core.supported._add;

pub fn asum(allocator: ?std.mem.Allocator, comptime T: type, x: NDArray(T)) !T {
    const supported = core.supported.whatSupportedNumericType(T);

    var res: T = undefined;
    switch (supported) {
        .BuiltinInt, .BuiltinFloat, .BuiltinBool => {
            res = 0;
        },
        .CustomComplexFloat => {
            res = T.init(0);
        },
        .CustomInt, .CustomReal, .CustomExpression => {
            if (!allocator) {
                return Error.NoAllocator;
            }
            res = T.init(allocator, 0);
        },
        .CustomComplex => {
            if (!allocator) {
                return Error.NoAllocator;
            }
            res = T.init(allocator, 0, 0);
        },
        .Unsupported => unreachable,
    }

    for (0..x.size) |i| {
        switch (supported) {
            .BuiltinInt, .BuiltinFloat, .BuiltinBool => {
                res += @abs(x.data[i]);
            },
            .CustomComplexFloat => {
                res.addInPlace(res, T.add(T.abs(x.data[i].re), T.abs(x.data[i].im)));
            },
            .CustomInt, .CustomReal, .CustomExpression => {
                res.addInPlace(T.abs(x.data[i]));
            },
            .CustomComplex => {
                res = 1; // Not done, all custom types to be reconsidered once they are implemented.
            },
            .Unsupported => unreachable,
        }
    }

    return res;
}

test "asum" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{ 8, 18, 7 }, .{});
    defer A.deinit();

    for (0..A.size) |i| {
        A.data[i] = 1;
    }

    const result = try NDArray(f64).BLAS.asum(null, A);
    try std.testing.expect(result == 1008);
}
