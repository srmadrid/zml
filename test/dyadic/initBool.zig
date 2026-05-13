const std = @import("std");

const zsl = @import("zsl");

test "initBool: true -> 1" {
    try std.testing.expectEqual(zsl.Dyadic(16, 16).one, zsl.Dyadic(16, 16).initValue(true));
}

test "initBool: false -> 0" {
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(false));
}

test "initBool: independent of mantissa/exponent widths" {
    inline for (.{
        .{ 1, 4 },  .{ 8, 4 },   .{ 16, 16 },
        .{ 32, 8 }, .{ 64, 32 }, .{ 128, 16 },
    }) |tc| {
        try std.testing.expectEqual(zsl.Dyadic(tc[0], tc[1]).one, zsl.Dyadic(tc[0], tc[1]).initValue(true));
        try std.testing.expectEqual(zsl.Dyadic(tc[0], tc[1]).zero, zsl.Dyadic(tc[0], tc[1]).initValue(false));
    }
}
