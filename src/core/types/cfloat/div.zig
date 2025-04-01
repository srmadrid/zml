const std = @import("std");
const types = @import("../../types.zig");
const Coerce = types.Coerce;

pub fn div(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (!types.isFixedPrecision(@TypeOf(left)) or !types.isFixedPrecision(@TypeOf(right)) or (!types.isComplex(@TypeOf(left)) and !types.isComplex(@TypeOf(right))))
        @compileError("At least one of left or right must be cfloat, the other must be an int, float, or cfloat");

    switch (types.numericType(@TypeOf(left))) {
        .int => {
            switch (types.numericType(@TypeOf(right))) {
                .int => unreachable,
                .float => unreachable,
                .cfloat => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}).inverse().mulReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), left, .{}));
                },
                else => unreachable,
            }
        },
        .float => {
            switch (types.numericType(@TypeOf(right))) {
                .int => unreachable,
                .float => unreachable,
                .cfloat => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}).inverse().mulReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), left, .{}));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(@TypeOf(right))) {
                .int => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}).divReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), right, .{}));
                },
                .float => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}).divReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), right, .{}));
                },
                .cfloat => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}).div(types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

test div {
    const x1: types.cf64 = .init(1, 2);
    const x2: types.cf64 = .init(3, 4);
    const x3: types.cf128 = .init(5, 6);
    const x4: f64 = 7;
    const x5: f128 = 8;

    const r1 = div(x1, x2);
    try std.testing.expectApproxEqRel(0.44, r1.re, 1e-16);
    try std.testing.expectApproxEqRel(0.08, r1.im, 1e-16);

    const r2 = div(x1, x3);
    try std.testing.expectApproxEqRel(0.278688524590163934426229508196721, r2.re, 1e-32);
    try std.testing.expectApproxEqRel(0.065573770491803278688524590163934, r2.im, 1e-32);

    const r3 = div(x1, x4);
    try std.testing.expectApproxEqRel(0.14285714285714285, r3.re, 1e-16);
    try std.testing.expectApproxEqRel(0.28571428571428571, r3.im, 1e-16);

    const r4 = div(x1, x5);
    try std.testing.expectApproxEqRel(0.125, r4.re, 1e-32);
    try std.testing.expectApproxEqRel(0.25, r4.im, 1e-32);

    const r5 = div(x3, x1);
    try std.testing.expectApproxEqRel(3.4, r5.re, 1e-32);
    try std.testing.expectApproxEqRel(-0.8, r5.im, 1e-32);

    const r6 = div(x3, x4);
    try std.testing.expectApproxEqRel(0.714285714285714285714285714285714, r6.re, 1e-32);
    try std.testing.expectApproxEqRel(0.857142857142857142857142857142857, r6.im, 1e-32);

    const r7 = div(x3, x5);
    try std.testing.expectApproxEqRel(0.625, r7.re, 1e-32);
    try std.testing.expectApproxEqRel(0.75, r7.im, 1e-32);

    const r8 = div(x4, x1);
    try std.testing.expectApproxEqRel(1.4, r8.re, 1e-16);
    try std.testing.expectApproxEqRel(-2.8, r8.im, 1e-16);

    const r9 = div(x4, x3);
    try std.testing.expectApproxEqRel(0.573770491803278688524590163934426, r9.re, 1e-32);
    try std.testing.expectApproxEqRel(-0.688524590163934426229508196721311, r9.im, 1e-32);

    const r10 = div(x5, x1);
    try std.testing.expectApproxEqRel(1.6, r10.re, 1e-32);
    try std.testing.expectApproxEqRel(-3.2, r10.im, 1e-32);

    const r11 = div(x5, x3);
    try std.testing.expectApproxEqRel(0.655737704918032786885245901639344, r11.re, 1e-32);
    try std.testing.expectApproxEqRel(-0.786885245901639344262295081967213, r11.im, 1e-32);
}
