const std = @import("std");
const zml = @import("zml");
const abs = zml.cfloat.abs;

const data_cf32: [25]struct { f32, zml.cf32 } = .{
    .{ 0xc.69ce3p+0, zml.cf32.init(0xcp-4, 0xc.64p+0) },
    .{ 0xc.69ce3p+0, zml.cf32.init(-0xc.64p+0, 0xcp-4) },
    .{ 0xc.69ce3p+0, zml.cf32.init(-0xcp-4, 0xc.64p+0) },
    .{ 0xc.69ce3p+0, zml.cf32.init(-0xc.64p+0, -0xcp-4) },
    .{ 0xc.69ce3p+0, zml.cf32.init(-0xcp-4, -0xc.64p+0) },
    .{ 0xcp-4, zml.cf32.init(-0xcp-4, 0x0p+0) },
    .{ 0xcp-4, zml.cf32.init(0xcp-4, 0x0p+0) },
    .{ 0x1p+0, zml.cf32.init(-0x1p+0, 0x0p+0) },
    .{ 0x1p+0, zml.cf32.init(0x1p+0, 0x0p+0) },
    .{ 0x3.65c04p+24, zml.cf32.init(-0x3.65c04p+24, 0x0p+0) },
    .{ 0x3.65c04p+24, zml.cf32.init(0x3.65c04p+24, 0x0p+0) },
    .{ 0x1.752e5p+0, zml.cf32.init(0xcp-4, 0x1.4p+0) },
    .{ 0xc.56714p+0, zml.cf32.init(-0x1.34be3p-4, -0xc.56623p+0) },
    .{ 0x1.2b0ffap+28, zml.cf32.init(-0x1.2b0ff8p+28, -0x2.549fc4p+16) },
    .{ 0x2.51109p-24, zml.cf32.init(-0x1.0932cp-80, -0x2.51109p-24) },
    .{ 0x1.055fb2p+48, zml.cf32.init(-0x1.055fb2p+48, 0x9.1ce86p+24) },
    .{ 0x1.26a566p+120, zml.cf32.init(-0x1.26a566p+120, 0x4.017b28p+92) },
    .{ 0x1.0eda54p+28, zml.cf32.init(-0x1.0eda54p+28, 0xb.09476p+0) },
    .{ 0xf.fffffp+124, zml.cf32.init(-0x1.133b84p+84, -0xf.fffffp+124) },
    .{ 0x3.4e5d78p+0, zml.cf32.init(-0x0p+0, -0x3.4e5d78p+0) },
    .{ 0x3.4e5d7cp+0, zml.cf32.init(-0x0p+0, -0x3.4e5d7cp+0) },
    .{ 0xa.21a95p+20, zml.cf32.init(-0xa.f59b8p+4, 0xa.21a95p+20) },
    .{ 0x1.e9d956p+56, zml.cf32.init(-0x1.30ed4cp+0, 0x1.e9d956p+56) },
    .{ 0x5.a5046p-4, zml.cf32.init(-0x1.250366p-36, -0x5.a5046p-4) },
    .{ 0x1.88858cp+84, zml.cf32.init(-0x1.88858cp+84, 0x5.bd9198p+36) },
};

const data_cf64: [27]struct { f64, zml.cf64 } = .{
    .{ 0xc.69ce375a71e08p+0, zml.cf64.init(0xcp-4, 0xc.64p+0) },
    .{ 0xc.69ce375a71e08p+0, zml.cf64.init(-0xc.64p+0, 0xcp-4) },
    .{ 0xc.69ce375a71e08p+0, zml.cf64.init(-0xcp-4, 0xc.64p+0) },
    .{ 0xc.69ce375a71e08p+0, zml.cf64.init(-0xc.64p+0, -0xcp-4) },
    .{ 0xc.69ce375a71e08p+0, zml.cf64.init(-0xcp-4, -0xc.64p+0) },
    .{ 0xcp-4, zml.cf64.init(-0xcp-4, 0x0p+0) },
    .{ 0xcp-4, zml.cf64.init(0xcp-4, 0x0p+0) },
    .{ 0x1p+0, zml.cf64.init(-0x1p+0, 0x0p+0) },
    .{ 0x1p+0, zml.cf64.init(0x1p+0, 0x0p+0) },
    .{ 0x3.65c04p+24, zml.cf64.init(-0x3.65c04p+24, 0x0p+0) },
    .{ 0x3.65c04p+24, zml.cf64.init(0x3.65c04p+24, 0x0p+0) },
    .{ 0x1.752e50db3a3a2p+0, zml.cf64.init(0xcp-4, 0x1.4p+0) },
    .{ 0xc.5671471794418p+0, zml.cf64.init(-0x1.34be3p-4, -0xc.56623p+0) },
    .{ 0x1.2b0ffa53208c7p+28, zml.cf64.init(-0x1.2b0ff8p+28, -0x2.549fc4p+16) },
    .{ 0x2.51109p-24, zml.cf64.init(-0x1.0932cp-80, -0x2.51109p-24) },
    .{ 0x1.055fb2000028bp+48, zml.cf64.init(-0x1.055fb2p+48, 0x9.1ce86p+24) },
    .{ 0x1.26a566p+120, zml.cf64.init(-0x1.26a566p+120, 0x4.017b28p+92) },
    .{ 0x1.0eda540000004p+28, zml.cf64.init(-0x1.0eda54p+28, 0xb.09476p+0) },
    .{ 0xf.fffffp+124, zml.cf64.init(-0x1.133b84p+84, -0xf.fffffp+124) },
    .{ 0xa.7d925f57f60cp+768, zml.cf64.init(-0x1.133b84p+84, -0xa.7d925f57f60cp+768) },
    .{ 0x3.4e5d78p+0, zml.cf64.init(-0x0p+0, -0x3.4e5d78p+0) },
    .{ 0x3.4e5d7cp+0, zml.cf64.init(-0x0p+0, -0x3.4e5d7cp+0) },
    .{ 0x3.4e5d7877324cp+0, zml.cf64.init(-0x0p+0, -0x3.4e5d7877324cp+0) },
    .{ 0xa.21a95005ed6f8p+20, zml.cf64.init(-0xa.f59b8p+4, 0xa.21a95p+20) },
    .{ 0x1.e9d956p+56, zml.cf64.init(-0x1.30ed4cp+0, 0x1.e9d956p+56) },
    .{ 0x5.a5046p-4, zml.cf64.init(-0x1.250366p-36, -0x5.a5046p-4) },
    .{ 0x1.88858cp+84, zml.cf64.init(-0x1.88858cp+84, 0x5.bd9198p+36) },
};

const data_cf80: [27]struct { f80, zml.cf80 } = .{
    .{ 0xc.69ce375a71e09aap+0, zml.cf80.init(0xcp-4, 0xc.64p+0) },
    .{ 0xc.69ce375a71e09aap+0, zml.cf80.init(-0xc.64p+0, 0xcp-4) },
    .{ 0xc.69ce375a71e09aap+0, zml.cf80.init(-0xcp-4, 0xc.64p+0) },
    .{ 0xc.69ce375a71e09aap+0, zml.cf80.init(-0xc.64p+0, -0xcp-4) },
    .{ 0xc.69ce375a71e09aap+0, zml.cf80.init(-0xcp-4, -0xc.64p+0) },
    .{ 0xcp-4, zml.cf80.init(-0xcp-4, 0x0p+0) },
    .{ 0xcp-4, zml.cf80.init(0xcp-4, 0x0p+0) },
    .{ 0x1p+0, zml.cf80.init(-0x1p+0, 0x0p+0) },
    .{ 0x1p+0, zml.cf80.init(0x1p+0, 0x0p+0) },
    .{ 0x3.65c04p+24, zml.cf80.init(-0x3.65c04p+24, 0x0p+0) },
    .{ 0x3.65c04p+24, zml.cf80.init(0x3.65c04p+24, 0x0p+0) },
    .{ 0x1.752e50db3a3a1b1cp+0, zml.cf80.init(0xcp-4, 0x1.4p+0) },
    .{ 0xc.56714717944142p+0, zml.cf80.init(-0x1.34be3p-4, -0xc.56623p+0) },
    .{ 0x1.2b0ffa53208c702cp+28, zml.cf80.init(-0x1.2b0ff8p+28, -0x2.549fc4p+16) },
    .{ 0x2.51109p-24, zml.cf80.init(-0x1.0932cp-80, -0x2.51109p-24) },
    .{ 0x1.055fb2000028ab42p+48, zml.cf80.init(-0x1.055fb2p+48, 0x9.1ce86p+24) },
    .{ 0x1.26a56600000006f8p+120, zml.cf80.init(-0x1.26a566p+120, 0x4.017b28p+92) },
    .{ 0x1.0eda54000000399p+28, zml.cf80.init(-0x1.0eda54p+28, 0xb.09476p+0) },
    .{ 0xf.fffffp+124, zml.cf80.init(-0x1.133b84p+84, -0xf.fffffp+124) },
    .{ 0xa.7d925f57f60cp+768, zml.cf80.init(-0x1.133b84p+84, -0xa.7d925f57f60cp+768) },
    .{ 0x3.4e5d78p+0, zml.cf80.init(-0x0p+0, -0x3.4e5d78p+0) },
    .{ 0x3.4e5d7cp+0, zml.cf80.init(-0x0p+0, -0x3.4e5d7cp+0) },
    .{ 0x3.4e5d7877324cp+0, zml.cf80.init(-0x0p+0, -0x3.4e5d7877324cp+0) },
    .{ 0xa.21a95005ed6fcp+20, zml.cf80.init(-0xa.f59b8p+4, 0xa.21a95p+20) },
    .{ 0x1.e9d956p+56, zml.cf80.init(-0x1.30ed4cp+0, 0x1.e9d956p+56) },
    .{ 0x5.a5046p-4, zml.cf80.init(-0x1.250366p-36, -0x5.a5046p-4) },
    .{ 0x1.88858cp+84, zml.cf80.init(-0x1.88858cp+84, 0x5.bd9198p+36) },
};

const data_cf128: [27]struct { f128, zml.cf128 } = .{
    .{ 0xc.69ce375a71e09a9df3616830c9e8p+0, zml.cf128.init(0xcp-4, 0xc.64p+0) },
    .{ 0xc.69ce375a71e09a9df3616830c9e8p+0, zml.cf128.init(-0xc.64p+0, 0xcp-4) },
    .{ 0xc.69ce375a71e09a9df3616830c9e8p+0, zml.cf128.init(-0xcp-4, 0xc.64p+0) },
    .{ 0xc.69ce375a71e09a9df3616830c9e8p+0, zml.cf128.init(-0xc.64p+0, -0xcp-4) },
    .{ 0xc.69ce375a71e09a9df3616830c9e8p+0, zml.cf128.init(-0xcp-4, -0xc.64p+0) },
    .{ 0xcp-4, zml.cf128.init(-0xcp-4, 0x0p+0) },
    .{ 0xcp-4, zml.cf128.init(0xcp-4, 0x0p+0) },
    .{ 0x1p+0, zml.cf128.init(-0x1p+0, 0x0p+0) },
    .{ 0x1p+0, zml.cf128.init(0x1p+0, 0x0p+0) },
    .{ 0x3.65c04p+24, zml.cf128.init(-0x3.65c04p+24, 0x0p+0) },
    .{ 0x3.65c04p+24, zml.cf128.init(0x3.65c04p+24, 0x0p+0) },
    .{ 0x1.752e50db3a3a1b1b33b0456f1fbbp+0, zml.cf128.init(0xcp-4, 0x1.4p+0) },
    .{ 0xc.56714717944141fc40fa4c791948p+0, zml.cf128.init(-0x1.34be3p-4, -0xc.56623p+0) },
    .{ 0x1.2b0ffa53208c702cbc8f252e1dfp+28, zml.cf128.init(-0x1.2b0ff8p+28, -0x2.549fc4p+16) },
    .{ 0x2.51109p-24, zml.cf128.init(-0x1.0932cp-80, -0x2.51109p-24) },
    .{ 0x1.055fb2000028ab411a37f7ed75bdp+48, zml.cf128.init(-0x1.055fb2p+48, 0x9.1ce86p+24) },
    .{ 0x1.26a56600000006f8887eefeb06d2p+120, zml.cf128.init(-0x1.26a566p+120, 0x4.017b28p+92) },
    .{ 0x1.0eda54000000398f4eef03909ac2p+28, zml.cf128.init(-0x1.0eda54p+28, 0xb.09476p+0) },
    .{ 0xf.fffff000000000000000093f4768p+124, zml.cf128.init(-0x1.133b84p+84, -0xf.fffffp+124) },
    .{ 0xa.7d925f57f60cp+768, zml.cf128.init(-0x1.133b84p+84, -0xa.7d925f57f60cp+768) },
    .{ 0x3.4e5d78p+0, zml.cf128.init(-0x0p+0, -0x3.4e5d78p+0) },
    .{ 0x3.4e5d7cp+0, zml.cf128.init(-0x0p+0, -0x3.4e5d7cp+0) },
    .{ 0x3.4e5d7877324cp+0, zml.cf128.init(-0x0p+0, -0x3.4e5d7877324cp+0) },
    .{ 0xa.21a95005ed6fbffe68d320c0fde8p+20, zml.cf128.init(-0xa.f59b8p+4, 0xa.21a95p+20) },
    .{ 0x1.e9d956p+56, zml.cf128.init(-0x1.30ed4cp+0, 0x1.e9d956p+56) },
    .{ 0x5.a5046000000000001db5376a4a08p-4, zml.cf128.init(-0x1.250366p-36, -0x5.a5046p-4) },
    .{ 0x1.88858c00000000000000000abef9p+84, zml.cf128.init(-0x1.88858cp+84, 0x5.bd9198p+36) },
};

test abs {
    for (data_cf32) |test_case| {
        try std.testing.expectEqual(test_case[0], abs(test_case[1]));
    }

    for (data_cf64) |test_case| {
        try std.testing.expectEqual(test_case[0], abs(test_case[1]));
    }

    for (data_cf80) |test_case| {
        try std.testing.expectEqual(test_case[0], abs(test_case[1]));
    }

    for (data_cf128) |test_case| {
        try std.testing.expectEqual(test_case[0], abs(test_case[1]));
    }
}
