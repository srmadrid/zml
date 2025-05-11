const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const cf128 = zml.cf128;
const sub = zml.cfloat.sub;

test sub {
    const x1: cf64 = .init(1, 2);
    const x2: cf64 = .init(3, 4);
    const x3: cf128 = .init(5, 6);
    const x4: f64 = 7;
    const x5: f128 = 8;

    const r1 = sub(x1, x2);
    try std.testing.expectEqual(cf64.init(-2, -2), r1);

    const r2 = sub(x1, x3);
    try std.testing.expectEqual(cf128.init(-4, -4), r2);

    const r3 = sub(x1, x4);
    try std.testing.expectEqual(cf64.init(-6, 2), r3);

    const r4 = sub(x1, x5);
    try std.testing.expectEqual(cf128.init(-7, 2), r4);

    const r5 = sub(x3, x1);
    try std.testing.expectEqual(cf128.init(4, 4), r5);

    const r6 = sub(x3, x4);
    try std.testing.expectEqual(cf128.init(-2, 6), r6);

    const r7 = sub(x3, x5);
    try std.testing.expectEqual(cf128.init(-3, 6), r7);

    const r8 = sub(x4, x1);
    try std.testing.expectEqual(cf64.init(6, -2), r8);

    const r9 = sub(x4, x3);
    try std.testing.expectEqual(cf128.init(2, -6), r9);

    const r10 = sub(x5, x1);
    try std.testing.expectEqual(cf128.init(7, -2), r10);

    const r11 = sub(x5, x3);
    try std.testing.expectEqual(cf128.init(3, -6), r11);
}
