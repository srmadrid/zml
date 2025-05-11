const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const cf128 = zml.cf128;
const mul = zml.cfloat.mul;

test mul {
    const x1: cf64 = .init(1, 2);
    const x2: cf64 = .init(3, 4);
    const x3: cf128 = .init(5, 6);
    const x4: f64 = 7;
    const x5: f128 = 8;

    const r1 = mul(x1, x2);
    try std.testing.expectEqual(cf64.init(-5, 10), r1);

    const r2 = mul(x1, x3);
    try std.testing.expectEqual(cf128.init(-7, 16), r2);

    const r3 = mul(x1, x4);
    try std.testing.expectEqual(cf64.init(7, 14), r3);

    const r4 = mul(x1, x5);
    try std.testing.expectEqual(cf128.init(8, 16), r4);

    const r5 = mul(x3, x1);
    try std.testing.expectEqual(cf128.init(-7, 16), r5);

    const r6 = mul(x3, x4);
    try std.testing.expectEqual(cf128.init(35, 42), r6);

    const r7 = mul(x3, x5);
    try std.testing.expectEqual(cf128.init(40, 48), r7);

    const r8 = mul(x4, x1);
    try std.testing.expectEqual(cf64.init(7, 14), r8);

    const r9 = mul(x4, x3);
    try std.testing.expectEqual(cf128.init(35, 42), r9);

    const r10 = mul(x5, x1);
    try std.testing.expectEqual(cf128.init(8, 16), r10);

    const r11 = mul(x5, x3);
    try std.testing.expectEqual(cf128.init(40, 48), r11);
}
