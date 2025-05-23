const std = @import("std");
const zml = @import("zml");
const cf32 = zml.cf32;
const cf64 = zml.cf64;
const cf80 = zml.cf80;
const cf128 = zml.cf128;
const arg = zml.cfloat.arg;

test arg {
    try std.testing.expectEqual(0x0p+0, arg(cf32.init(std.math.inf(f32), 0x2p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf32.init(std.math.inf(f32), -0x2p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, arg(cf32.init(0xap+0, std.math.inf(f32))));
    try std.testing.expectEqual(-0x1.921fb6p+0, arg(cf32.init(0xap+0, -std.math.inf(f32))));
    try std.testing.expectEqual(0x3.243f6cp+0, arg(cf32.init(-std.math.inf(f32), 0xap+0)));
    try std.testing.expectEqual(-0x3.243f6cp+0, arg(cf32.init(-std.math.inf(f32), -0xap+0)));
    try std.testing.expectEqual(0xc.90fdbp-4, arg(cf32.init(std.math.inf(f32), std.math.inf(f32))));
    try std.testing.expectEqual(-0xc.90fdbp-4, arg(cf32.init(std.math.inf(f32), -std.math.inf(f32))));
    try std.testing.expectEqual(0x2.5b2f9p+0, arg(cf32.init(-std.math.inf(f32), std.math.inf(f32))));
    try std.testing.expectEqual(-0x2.5b2f9p+0, arg(cf32.init(-std.math.inf(f32), -std.math.inf(f32))));
    try std.testing.expectEqual(0x0p+0, arg(cf32.init(0x2p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf32.init(0x2p+0, -0x0p+0)));
    try std.testing.expectEqual(0x0p+0, arg(cf32.init(0x0p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf32.init(0x0p+0, -0x0p+0)));
    try std.testing.expectEqual(0x3.243f6cp+0, arg(cf32.init(-0x2p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6cp+0, arg(cf32.init(-0x2p+0, -0x0p+0)));
    try std.testing.expectEqual(0x3.243f6cp+0, arg(cf32.init(-0x0p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6cp+0, arg(cf32.init(-0x0p+0, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, arg(cf32.init(0x0p+0, 0x2p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, arg(cf32.init(-0x0p+0, 0x2p+0)));
    try std.testing.expectEqual(-0x1.921fb6p+0, arg(cf32.init(0x0p+0, -0x2p+0)));
    try std.testing.expectEqual(-0x1.921fb6p+0, arg(cf32.init(-0x0p+0, -0x2p+0)));
    try std.testing.expectEqual(0x1.9c22cep-4, arg(cf32.init(0x2.f2f308p+0, 0x4.c3841p-4)));
    try std.testing.expectEqual(-0x1.1dd4c4p-4, arg(cf32.init(0xd.3de7ap-36, -0xe.cf143p-40)));
    try std.testing.expectEqual(0x2.7c1784p-4, arg(cf32.init(0x2.21e65p+0, 0x5.576cf8p-4)));
    try std.testing.expectEqual(-0x2.1dbac4p-4, arg(cf32.init(0x1.f4755cp+0, -0x4.29411p-4)));
    try std.testing.expectEqual(-0x1.921fb6p+0, arg(cf32.init(-0xf.9c4c8p-4, -0xa.b4101p+20)));
    try std.testing.expectEqual(0x9.23e97p-8, arg(cf32.init(0x7.40ac68p+0, 0x4.251bb8p-4)));
    try std.testing.expectEqual(0x1.fd0a44p-4, arg(cf32.init(0xa.3ac3cp+68, 0x1.47239ep+68)));
    try std.testing.expectEqual(-0x1.de8936p-4, arg(cf32.init(0x3.8ff10cp+0, -0x6.b0794p-4)));
    try std.testing.expectEqual(-0x1.921fb6p+0, arg(cf32.init(-0x3.973cc4p+72, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.9220ccp+0, arg(cf32.init(-0x1.0a512ap-120, 0xf.54681p-108)));
    try std.testing.expectEqual(-0x1.e5bac4p+0, arg(cf32.init(-0x3.be0054p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(-0x1.e5bac4p+0, arg(cf32.init(-0x3.be0058p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(0x1.921fb6p+0, arg(cf32.init(-0x1.0236b6p-20, 0x2.a6e504p+108)));
    try std.testing.expectEqual(0x3.2ec6ep-4, arg(cf32.init(0x9.27b6p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6dcp-4, arg(cf32.init(0x9.27b6p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e8p-4, arg(cf32.init(0x9.27b5fp+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6e4p-4, arg(cf32.init(0x9.27b5fp+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0xc.90fdbp-4, arg(cf32.init(0x8p-152, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, arg(cf32.init(0x8p-152, 0x0p+0)));

    try std.testing.expectEqual(0x0p+0, arg(cf64.init(std.math.inf(f64), 0x2p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf64.init(std.math.inf(f64), -0x2p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, arg(cf64.init(0xap+0, std.math.inf(f64))));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, arg(cf64.init(0xap+0, -std.math.inf(f64))));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, arg(cf64.init(-std.math.inf(f64), 0xap+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, arg(cf64.init(-std.math.inf(f64), -0xap+0)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, arg(cf64.init(std.math.inf(f64), std.math.inf(f64))));
    try std.testing.expectEqual(-0xc.90fdaa22168cp-4, arg(cf64.init(std.math.inf(f64), -std.math.inf(f64))));
    try std.testing.expectEqual(0x2.5b2f8fe6643a4p+0, arg(cf64.init(-std.math.inf(f64), std.math.inf(f64))));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a4p+0, arg(cf64.init(-std.math.inf(f64), -std.math.inf(f64))));
    try std.testing.expectEqual(0x0p+0, arg(cf64.init(0x2p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf64.init(0x2p+0, -0x0p+0)));
    try std.testing.expectEqual(0x0p+0, arg(cf64.init(0x0p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf64.init(0x0p+0, -0x0p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, arg(cf64.init(-0x2p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, arg(cf64.init(-0x2p+0, -0x0p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, arg(cf64.init(-0x0p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, arg(cf64.init(-0x0p+0, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, arg(cf64.init(0x0p+0, 0x2p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, arg(cf64.init(-0x0p+0, 0x2p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, arg(cf64.init(0x0p+0, -0x2p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, arg(cf64.init(-0x0p+0, -0x2p+0)));
    try std.testing.expectEqual(0x1.9c22ce44a722ap-4, arg(cf64.init(0x2.f2f308p+0, 0x4.c3841p-4)));
    try std.testing.expectEqual(-0x1.1dd4c4e264577p-4, arg(cf64.init(0xd.3de7ap-36, -0xe.cf143p-40)));
    try std.testing.expectEqual(0x2.7c1782a75e16cp-4, arg(cf64.init(0x2.21e65p+0, 0x5.576cf8p-4)));
    try std.testing.expectEqual(-0x2.1dbac4fa4bfecp-4, arg(cf64.init(0x1.f4755cp+0, -0x4.29411p-4)));
    try std.testing.expectEqual(-0x1.921fb6b9a118dp+0, arg(cf64.init(-0xf.9c4c8p-4, -0xa.b4101p+20)));
    try std.testing.expectEqual(0x9.23e97736442d8p-8, arg(cf64.init(0x7.40ac68p+0, 0x4.251bb8p-4)));
    try std.testing.expectEqual(0x1.fd0a44d0aba44p-4, arg(cf64.init(0xa.3ac3cp+68, 0x1.47239ep+68)));
    try std.testing.expectEqual(-0x1.de89352a0e839p-4, arg(cf64.init(0x3.8ff10cp+0, -0x6.b0794p-4)));
    try std.testing.expectEqual(-0x1.921fb54442d19p+0, arg(cf64.init(-0x3.973cc4p+72, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.9220cb3a73868p+0, arg(cf64.init(-0x1.0a512ap-120, 0xf.54681p-108)));
    try std.testing.expectEqual(-0x1.e5bac45eb390ap+0, arg(cf64.init(-0x3.be0054p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(-0x1.e5bac4b1d8c3dp+0, arg(cf64.init(-0x3.be0058p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(-0x1.e5bac46572919p+0, arg(cf64.init(-0x3.be0054531569p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, arg(cf64.init(-0x1.0236b6p-20, 0x2.a6e504p+108)));
    try std.testing.expectEqual(0x3.2ec6e0e7c264p-4, arg(cf64.init(0x9.27b6p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6dd8be6cc2p-4, arg(cf64.init(0x9.27b6p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e02be7d16p-4, arg(cf64.init(0x9.27b6p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e652712d4p-4, arg(cf64.init(0x9.27b5fp+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6e2f6959p-4, arg(cf64.init(0x9.27b5fp+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e59696998p-4, arg(cf64.init(0x9.27b5fp+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea68p-4, arg(cf64.init(0x9.27b5fd9157b7p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30dep-4, arg(cf64.init(0x9.27b5fd9157b7p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e0feb413cp-4, arg(cf64.init(0x9.27b5fd9157b7p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea6cp-4, arg(cf64.init(0x9.27b5fd9157b68p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30ep-4, arg(cf64.init(0x9.27b5fd9157b68p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e0feb413ep-4, arg(cf64.init(0x9.27b5fd9157b68p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, arg(cf64.init(0x8p-152, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, arg(cf64.init(0x8p-152, 0x0p+0)));
    try std.testing.expectEqual(0x8p-928, arg(cf64.init(0x8p-152, 0x4p-1076)));

    try std.testing.expectEqual(0x0p+0, arg(cf80.init(std.math.inf(f80), 0x2p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf80.init(std.math.inf(f80), -0x2p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, arg(cf80.init(0xap+0, std.math.inf(f80))));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, arg(cf80.init(0xap+0, -std.math.inf(f80))));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, arg(cf80.init(-std.math.inf(f80), 0xap+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, arg(cf80.init(-std.math.inf(f80), -0xap+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, arg(cf80.init(std.math.inf(f80), std.math.inf(f80))));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, arg(cf80.init(std.math.inf(f80), -std.math.inf(f80))));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, arg(cf80.init(-std.math.inf(f80), std.math.inf(f80))));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, arg(cf80.init(-std.math.inf(f80), -std.math.inf(f80))));
    try std.testing.expectEqual(0x0p+0, arg(cf80.init(0x2p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf80.init(0x2p+0, -0x0p+0)));
    try std.testing.expectEqual(0x0p+0, arg(cf80.init(0x0p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf80.init(0x0p+0, -0x0p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, arg(cf80.init(-0x2p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, arg(cf80.init(-0x2p+0, -0x0p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, arg(cf80.init(-0x0p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, arg(cf80.init(-0x0p+0, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, arg(cf80.init(0x0p+0, 0x2p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, arg(cf80.init(-0x0p+0, 0x2p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, arg(cf80.init(0x0p+0, -0x2p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, arg(cf80.init(-0x0p+0, -0x2p+0)));
    try std.testing.expectEqual(0x1.9c22ce44a7229d12p-4, arg(cf80.init(0x2.f2f308p+0, 0x4.c3841p-4)));
    try std.testing.expectEqual(-0x1.1dd4c4e2645769d2p-4, arg(cf80.init(0xd.3de7ap-36, -0xe.cf143p-40)));
    try std.testing.expectEqual(0x2.7c1782a75e16b744p-4, arg(cf80.init(0x2.21e65p+0, 0x5.576cf8p-4)));
    try std.testing.expectEqual(-0x2.1dbac4fa4bfeb75p-4, arg(cf80.init(0x1.f4755cp+0, -0x4.29411p-4)));
    try std.testing.expectEqual(-0x1.921fb6b9a118c896p+0, arg(cf80.init(-0xf.9c4c8p-4, -0xa.b4101p+20)));
    try std.testing.expectEqual(0x9.23e97736442d916p-8, arg(cf80.init(0x7.40ac68p+0, 0x4.251bb8p-4)));
    try std.testing.expectEqual(0x1.fd0a44d0aba440f4p-4, arg(cf80.init(0xa.3ac3cp+68, 0x1.47239ep+68)));
    try std.testing.expectEqual(-0x1.de89352a0e839634p-4, arg(cf80.init(0x3.8ff10cp+0, -0x6.b0794p-4)));
    try std.testing.expectEqual(-0x1.921fb54442d188p+0, arg(cf80.init(-0x3.973cc4p+72, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.9220cb3a738682a6p+0, arg(cf80.init(-0x1.0a512ap-120, 0xf.54681p-108)));
    try std.testing.expectEqual(-0x1.e5bac45eb390a6d4p+0, arg(cf80.init(-0x3.be0054p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(-0x1.e5bac4b1d8c3c81ap+0, arg(cf80.init(-0x3.be0058p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(-0x1.e5bac46572919652p+0, arg(cf80.init(-0x3.be0054531569p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, arg(cf80.init(-0x1.0236b6p-20, 0x2.a6e504p+108)));
    try std.testing.expectEqual(0x3.2ec6e0e7c2640134p-4, arg(cf80.init(0x9.27b6p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6dd8be6cc206p-4, arg(cf80.init(0x9.27b6p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e02be7d167c8p-4, arg(cf80.init(0x9.27b6p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e652712d3ebcp-4, arg(cf80.init(0x9.27b5fp+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6e2f6958ff48cp-4, arg(cf80.init(0x9.27b5fp+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e596969976a4p-4, arg(cf80.init(0x9.27b5fp+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea68638p-4, arg(cf80.init(0x9.27b5fd9157b7p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30dd2ccp-4, arg(cf80.init(0x9.27b5fd9157b7p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e0feb413bec8p-4, arg(cf80.init(0x9.27b5fd9157b7p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea6b18cp-4, arg(cf80.init(0x9.27b5fd9157b68p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30dfe2p-4, arg(cf80.init(0x9.27b5fd9157b68p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e0feb413ea1cp-4, arg(cf80.init(0x9.27b5fd9157b68p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea698c4p-4, arg(cf80.init(0x9.27b5fd9157b6c93p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30de558p-4, arg(cf80.init(0x9.27b5fd9157b6c93p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e0feb413d154p-4, arg(cf80.init(0x9.27b5fd9157b6c93p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea698ccp-4, arg(cf80.init(0x9.27b5fd9157b6c92p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30de56p-4, arg(cf80.init(0x9.27b5fd9157b6c92p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e0feb413d15cp-4, arg(cf80.init(0x9.27b5fd9157b6c92p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, arg(cf80.init(0x8p-152, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, arg(cf80.init(0x8p-152, 0x0p+0)));
    try std.testing.expectEqual(0x8p-928, arg(cf80.init(0x8p-152, 0x4p-1076)));

    try std.testing.expectEqual(0x0p+0, arg(cf128.init(std.math.inf(f128), 0x2p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf128.init(std.math.inf(f128), -0x2p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, arg(cf128.init(0xap+0, std.math.inf(f128))));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, arg(cf128.init(0xap+0, -std.math.inf(f128))));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, arg(cf128.init(-std.math.inf(f128), 0xap+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, arg(cf128.init(-std.math.inf(f128), -0xap+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, arg(cf128.init(std.math.inf(f128), std.math.inf(f128))));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, arg(cf128.init(std.math.inf(f128), -std.math.inf(f128))));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, arg(cf128.init(-std.math.inf(f128), std.math.inf(f128))));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, arg(cf128.init(-std.math.inf(f128), -std.math.inf(f128))));
    try std.testing.expectEqual(0x0p+0, arg(cf128.init(0x2p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf128.init(0x2p+0, -0x0p+0)));
    try std.testing.expectEqual(0x0p+0, arg(cf128.init(0x0p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, arg(cf128.init(0x0p+0, -0x0p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, arg(cf128.init(-0x2p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, arg(cf128.init(-0x2p+0, -0x0p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, arg(cf128.init(-0x0p+0, 0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, arg(cf128.init(-0x0p+0, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, arg(cf128.init(0x0p+0, 0x2p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, arg(cf128.init(-0x0p+0, 0x2p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, arg(cf128.init(0x0p+0, -0x2p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, arg(cf128.init(-0x0p+0, -0x2p+0)));
    // try std.testing.expectEqual(0x1.9c22ce44a7229d114c2b882266fap-4, arg(cf128.init(0x2.f2f308p+0, 0x4.c3841p-4)));
    try std.testing.expectEqual(-0x1.1dd4c4e2645769d1f7ebdc32a451p-4, arg(cf128.init(0xd.3de7ap-36, -0xe.cf143p-40)));
    try std.testing.expectEqual(0x2.7c1782a75e16b743e48c247c62cap-4, arg(cf128.init(0x2.21e65p+0, 0x5.576cf8p-4)));
    try std.testing.expectEqual(-0x2.1dbac4fa4bfeb74f6140009955a6p-4, arg(cf128.init(0x1.f4755cp+0, -0x4.29411p-4)));
    try std.testing.expectEqual(-0x1.921fb6b9a118c89590d474178551p+0, arg(cf128.init(-0xf.9c4c8p-4, -0xa.b4101p+20)));
    // try std.testing.expectEqual(0x9.23e97736442d915917b21858b148p-8, arg(cf128.init(0x7.40ac68p+0, 0x4.251bb8p-4)));
    // try std.testing.expectEqual(0x1.fd0a44d0aba440f30193e8545bc2p-4, arg(cf128.init(0xa.3ac3cp+68, 0x1.47239ep+68)));
    try std.testing.expectEqual(-0x1.de89352a0e839633c32d65e25422p-4, arg(cf128.init(0x3.8ff10cp+0, -0x6.b0794p-4)));
    try std.testing.expectEqual(-0x1.921fb54442d18800c6545c53c94fp+0, arg(cf128.init(-0x3.973cc4p+72, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.9220cb3a738682a53ab7ff520774p+0, arg(cf128.init(-0x1.0a512ap-120, 0xf.54681p-108)));
    try std.testing.expectEqual(-0x1.e5bac45eb390a6d33f541e8704aap+0, arg(cf128.init(-0x3.be0054p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(-0x1.e5bac4b1d8c3c81987f3ee7bd019p+0, arg(cf128.init(-0x3.be0058p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(-0x1.e5bac4657291965130fd9b80bae4p+0, arg(cf128.init(-0x3.be0054531569p-4, -0xb.0c5a9p-4)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, arg(cf128.init(-0x1.0236b6p-20, 0x2.a6e504p+108)));
    try std.testing.expectEqual(0x3.2ec6e0e7c2640133b01126200202p-4, arg(cf128.init(0x9.27b6p+0, 0x1.d8759cp+0)));
    // try std.testing.expectEqual(0x3.2ec6dd8be6cc20610e9a6ce98726p-4, arg(cf128.init(0x9.27b6p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e02be7d167c686712d2dadd8p-4, arg(cf128.init(0x9.27b6p+0, 0x1.d8759b9024992p+0)));
    // try std.testing.expectEqual(0x3.2ec6e652712d3ebc292bc227adf8p-4, arg(cf128.init(0x9.27b5fp+0, 0x1.d8759cp+0)));
    // try std.testing.expectEqual(0x3.2ec6e2f6958ff48a5bba9903412ap-4, arg(cf128.init(0x9.27b5fp+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e596969976a535107fa856ecp-4, arg(cf128.init(0x9.27b5fp+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea6863714534ed3f2c4p-4, arg(cf128.init(0x9.27b5fd9157b7p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30dd2cb350960010856p-4, arg(cf128.init(0x9.27b5fd9157b7p+0, 0x1.d8759ap+0)));
    // try std.testing.expectEqual(0x3.2ec6e0feb413bec772f42133ed3ap-4, arg(cf128.init(0x9.27b5fd9157b7p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea6b18c8a6a94e18e8ep-4, arg(cf128.init(0x9.27b5fd9157b68p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30dfe20aaf55b15729cp-4, arg(cf128.init(0x9.27b5fd9157b68p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e0feb413ea1ce901f1f33f48p-4, arg(cf128.init(0x9.27b5fd9157b68p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea698c58d050633860cp-4, arg(cf128.init(0x9.27b5fd9157b6c93p+0, 0x1.d8759cp+0)));
    // try std.testing.expectEqual(0x3.2ec6de5eb30de559ada88d65a5cap-4, arg(cf128.init(0x9.27b5fd9157b6c93p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e0feb413d155eba1cbb7da3p-4, arg(cf128.init(0x9.27b5fd9157b6c93p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea698caf7b3c91c47cp-4, arg(cf128.init(0x9.27b5fd9157b6c92p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30de55f18574ae50858p-4, arg(cf128.init(0x9.27b5fd9157b6c92p+0, 0x1.d8759ap+0)));
    // try std.testing.expectEqual(0x3.2ec6e0feb413d15b56508d71f21ap-4, arg(cf128.init(0x9.27b5fd9157b6c92p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea698c737390887086p-4, arg(cf128.init(0x9.27b5fd9157b6c92b151371ca23d8p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30de55b57dc8e0f5b42p-4, arg(cf128.init(0x9.27b5fd9157b6c92b151371ca23d8p+0, 0x1.d8759ap+0)));
    // try std.testing.expectEqual(0x3.2ec6e0feb413d15795d5cdae5624p-4, arg(cf128.init(0x9.27b5fd9157b6c92b151371ca23d8p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0x3.2ec6e1ba8ea698c7373908870852p-4, arg(cf128.init(0x9.27b5fd9157b6c92b151371ca24p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30de55b57dc8e0f5b34p-4, arg(cf128.init(0x9.27b5fd9157b6c92b151371ca24p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e0feb413d15795d5cdae5616p-4, arg(cf128.init(0x9.27b5fd9157b6c92b151371ca24p+0, 0x1.d8759b9024992p+0)));
    // try std.testing.expectEqual(0x3.2ec6e1ba8ea698c73739088709aep-4, arg(cf128.init(0x9.27b5fd9157b6c92b151371ca2p+0, 0x1.d8759cp+0)));
    try std.testing.expectEqual(0x3.2ec6de5eb30de55b57dc8e0f5c9p-4, arg(cf128.init(0x9.27b5fd9157b6c92b151371ca2p+0, 0x1.d8759ap+0)));
    try std.testing.expectEqual(0x3.2ec6e0feb413d15795d5cdae577p-4, arg(cf128.init(0x9.27b5fd9157b6c92b151371ca2p+0, 0x1.d8759b9024992p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, arg(cf128.init(0x8p-152, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, arg(cf128.init(0x8p-152, 0x0p+0)));
    try std.testing.expectEqual(0x8p-928, arg(cf128.init(0x8p-152, 0x4p-1076)));
}
