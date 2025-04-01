const std = @import("std");
const types = @import("../../types.zig");
const Coerce = types.Coerce;

pub fn add(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (!types.isFixedPrecision(@TypeOf(left)) or !types.isFixedPrecision(@TypeOf(right)) or (!types.isComplex(@TypeOf(left)) and !types.isComplex(@TypeOf(right))))
        @compileError("At least one of left or right must be cfloat, the other must be an int, float, or cfloat");

    switch (types.numericType(@TypeOf(left))) {
        .int => {
            switch (types.numericType(@TypeOf(right))) {
                .int => unreachable,
                .float => unreachable,
                .cfloat => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}).addReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), left, .{}));
                },
                else => unreachable,
            }
        },
        .float => {
            switch (types.numericType(@TypeOf(right))) {
                .int => unreachable,
                .float => unreachable,
                .cfloat => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}).addReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), left, .{}));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(@TypeOf(right))) {
                .int => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}).addReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), right, .{}));
                },
                .float => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}).addReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), right, .{}));
                },
                .cfloat => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}).add(types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

test add {
    const x1: types.cf64 = .init(1, 2);
    const x2: types.cf64 = .init(3, 4);
    const x3: types.cf128 = .init(5, 6);
    const x4: f64 = 7;
    const x5: f128 = 8;

    const r1 = add(x1, x2);
    try std.testing.expectEqual(types.cf64.init(4, 6), r1);

    const r2 = add(x1, x3);
    try std.testing.expectEqual(types.cf128.init(6, 8), r2);

    const r3 = add(x1, x4);
    try std.testing.expectEqual(types.cf64.init(8, 2), r3);

    const r4 = add(x1, x5);
    try std.testing.expectEqual(types.cf128.init(9, 2), r4);

    const r5 = add(x3, x1);
    try std.testing.expectEqual(types.cf128.init(6, 8), r5);

    const r6 = add(x3, x4);
    try std.testing.expectEqual(types.cf128.init(12, 6), r6);

    const r7 = add(x3, x5);
    try std.testing.expectEqual(types.cf128.init(13, 6), r7);

    const r8 = add(x4, x1);
    try std.testing.expectEqual(types.cf64.init(8, 2), r8);

    const r9 = add(x4, x3);
    try std.testing.expectEqual(types.cf128.init(12, 6), r9);

    const r10 = add(x5, x1);
    try std.testing.expectEqual(types.cf128.init(9, 2), r10);

    const r11 = add(x5, x3);
    try std.testing.expectEqual(types.cf128.init(13, 6), r11);
}
