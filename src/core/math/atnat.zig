const endian = @import("builtin").cpu.arch.endian();

pub const HIGH_HALF: u32 = if (endian == .big) 0 else 1;
pub const LOW_HALF: u32 = if (endian == .big) 1 else 0;
pub const d3: [2]u32 = if (endian == .big) .{ 0xbfd55555, 0x55555555 } else .{ 0x55555555, 0xbfd55555 }; // -0.333...
pub const d5: [2]u32 = if (endian == .big) .{ 0x3fc99999, 0x999997fd } else .{ 0x999997fd, 0x3fc99999 }; // 0.199...
pub const d7: [2]u32 = if (endian == .big) .{ 0xbfc24924, 0x923f7603 } else .{ 0x923f7603, 0xbfc24924 }; // -0.142...
pub const d9: [2]u32 = if (endian == .big) .{ 0x3fbc71c6, 0xe5129a3b } else .{ 0xe5129a3b, 0x3fbc71c6 }; // 0.111...
pub const d11: [2]u32 = if (endian == .big) .{ 0xbfb74580, 0x22b13c25 } else .{ 0x22b13c25, 0xbfb74580 }; // -0.090...
pub const d13: [2]u32 = if (endian == .big) .{ 0x3fb375f0, 0x8b31cbce } else .{ 0x8b31cbce, 0x3fb375f0 }; // 0.076...
pub const f3: [2]u32 = if (endian == .big) .{ 0xbfd55555, 0x55555555 } else .{ 0x55555555, 0xbfd55555 }; // -1/3
pub const ff3: [2]u32 = if (endian == .big) .{ 0xbc755555, 0x55555555 } else .{ 0x55555555, 0xbc755555 }; // -1/3-f3
pub const f5: [2]u32 = if (endian == .big) .{ 0x3fc99999, 0x9999999a } else .{ 0x9999999a, 0x3fc99999 }; // 1/5
pub const ff5: [2]u32 = if (endian == .big) .{ 0xbc699999, 0x9999999a } else .{ 0x9999999a, 0xbc699999 }; // 1/5-f5
pub const f7: [2]u32 = if (endian == .big) .{ 0xbfc24924, 0x92492492 } else .{ 0x92492492, 0xbfc24924 }; // -1/7
pub const ff7: [2]u32 = if (endian == .big) .{ 0xbc624924, 0x92492492 } else .{ 0x92492492, 0xbc624924 }; // -1/7-f7
pub const f9: [2]u32 = if (endian == .big) .{ 0x3fbc71c7, 0x1c71c71c } else .{ 0x1c71c71c, 0x3fbc71c7 }; // 1/9
pub const ff9: [2]u32 = if (endian == .big) .{ 0x3c5c71c7, 0x1c71c71c } else .{ 0x1c71c71c, 0x3c5c71c7 }; // 1/9-f9
pub const f11: [2]u32 = if (endian == .big) .{ 0xbfb745d1, 0x745d1746 } else .{ 0x745d1746, 0xbfb745d1 }; // -1/11
pub const f13: [2]u32 = if (endian == .big) .{ 0x3fb3b13b, 0x13b13b14 } else .{ 0x13b13b14, 0x3fb3b13b }; // 1/13
pub const f15: [2]u32 = if (endian == .big) .{ 0xbfb11111, 0x11111111 } else .{ 0x11111111, 0xbfb11111 }; // -1/15
pub const f17: [2]u32 = if (endian == .big) .{ 0x3fae1e1e, 0x1e1e1e1e } else .{ 0x1e1e1e1e, 0x3fae1e1e }; // 1/17
pub const f19: [2]u32 = if (endian == .big) .{ 0xbfaaf286, 0xbca1af28 } else .{ 0xbca1af28, 0xbfaaf286 }; // -1/19
pub const a: [2]u32 = if (endian == .big) .{ 0x3e4bb67a, 0x00000000 } else .{ 0x00000000, 0x3e4bb67a }; // 1.290e-8
pub const b: [2]u32 = if (endian == .big) .{ 0x3fb00000, 0x00000000 } else .{ 0x00000000, 0x3fb00000 }; // 1/16
pub const c: [2]u32 = if (endian == .big) .{ 0x3ff00000, 0x00000000 } else .{ 0x00000000, 0x3ff00000 }; // 1
pub const d: [2]u32 = if (endian == .big) .{ 0x40300000, 0x00000000 } else .{ 0x00000000, 0x40300000 }; // 16
pub const e: [2]u32 = if (endian == .big) .{ 0x43349ff2, 0x00000000 } else .{ 0x00000000, 0x43349ff2 }; // 5.805e15
pub const hpi: [2]u32 = if (endian == .big) .{ 0x3ff921fb, 0x54442d18 } else .{ 0x54442d18, 0x3ff921fb }; // pi/2
pub const mhpi: [2]u32 = if (endian == .big) .{ 0xbff921fb, 0x54442d18 } else .{ 0x54442d18, 0xbff921fb }; // -pi/2
pub const hpi1: [2]u32 = if (endian == .big) .{ 0x3c91a626, 0x33145c07 } else .{ 0x33145c07, 0x3c91a626 }; // pi/2-hpi
pub const U1: [2]u32 = if (endian == .big) .{ 0x3c2d3382, 0x00000000 } else .{ 0x00000000, 0x3c2d3382 }; // 7.915e-19
pub const U21: [2]u32 = if (endian == .big) .{ 0x3c6dffc0, 0x00000000 } else .{ 0x00000000, 0x3c6dffc0 }; // 1.301e-17
pub const U22: [2]u32 = if (endian == .big) .{ 0x3c527bd0, 0x00000000 } else .{ 0x00000000, 0x3c527bd0 }; // 4.008e-18
pub const U23: [2]u32 = if (endian == .big) .{ 0x3c3cd057, 0x00000000 } else .{ 0x00000000, 0x3c3cd057 }; // 1.562e-18
pub const U24: [2]u32 = if (endian == .big) .{ 0x3c329cdf, 0x00000000 } else .{ 0x00000000, 0x3c329cdf }; // 1.009e-18
pub const U31: [2]u32 = if (endian == .big) .{ 0x3c3a1edf, 0x00000000 } else .{ 0x00000000, 0x3c3a1edf }; // 1.416e-18
pub const U32: [2]u32 = if (endian == .big) .{ 0x3c33f0e1, 0x00000000 } else .{ 0x00000000, 0x3c33f0e1 }; // 1.081e-18
pub const U4: [2]u32 = if (endian == .big) .{ 0x3bf955e4, 0x00000000 } else .{ 0x00000000, 0x3bf955e4 }; // 8.584e-20
pub const U5: [2]u32 = if (endian == .big) .{ 0x3aaef2d1, 0x00000000 } else .{ 0x00000000, 0x3aaef2d1 }; // 5e-26
pub const U6: [2]u32 = if (endian == .big) .{ 0x3a98c56d, 0x00000000 } else .{ 0x00000000, 0x3a98c56d }; // 2.001e-26
pub const U7: [2]u32 = if (endian == .big) .{ 0x3a9375de, 0x00000000 } else .{ 0x00000000, 0x3a9375de }; // 1.572e-26
pub const U8: [2]u32 = if (endian == .big) .{ 0x3a6eeb36, 0x00000000 } else .{ 0x00000000, 0x3a6eeb36 }; // 3.122e-27
pub const U9: [4][2]u32 = if (endian == .big) .{ // 2.658e-35
    .{ 0x38c1aa5b, 0x00000000 }, // 2.658e-35
    .{ 0x35c1aa4d, 0x00000000 }, // 9.443e-50
    .{ 0x32c1aa88, 0x00000000 }, // 3.355e-64
    .{ 0x11c1aa56, 0x00000000 }, // 3.818e-223
} else .{
    .{ 0x00000000, 0x38c1aa5b }, // 2.658e-35
    .{ 0x00000000, 0x35c1aa4d }, // 9.443e-50
    .{ 0x00000000, 0x32c1aa88 }, // 3.355e-64
    .{ 0x00000000, 0x11c1aa56 }, // 3.818e-223
};

pub const cij = @import("uatan_tbl.zig").cij;
