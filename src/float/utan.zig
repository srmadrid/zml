const endian = @import("builtin").cpu.arch.endian();

pub const HIGH_HALF: u32 = if (endian == .big) 0 else 1;
pub const LOW_HALF: u32 = if (endian == .big) 1 else 0;
// polynomial I
pub const d3: [2]u32 = if (endian == .big) .{ 0x3fd55555, 0x55555555 } else .{ 0x55555555, 0x3fd55555 }; //  0.333...
pub const d5: [2]u32 = if (endian == .big) .{ 0x3fc11111, 0x111107c6 } else .{ 0x111107c6, 0x3fc11111 }; //  0.133...
pub const d7: [2]u32 = if (endian == .big) .{ 0x3faba1ba, 0x1cdb8745 } else .{ 0x1cdb8745, 0x3faba1ba }; //    .
pub const d9: [2]u32 = if (endian == .big) .{ 0x3f9664ed, 0x49cfc666 } else .{ 0x49cfc666, 0x3f9664ed }; //    .
pub const d11: [2]u32 = if (endian == .big) .{ 0x3f82385a, 0x3cf2e4ea } else .{ 0x3cf2e4ea, 0x3f82385a }; //    .
// polynomial II
// polynomial III
pub const e0: [2]u32 = if (endian == .big) .{ 0x3fd55555, 0x55554dbd } else .{ 0x55554dbd, 0x3fd55555 }; //    .
pub const e1: [2]u32 = if (endian == .big) .{ 0x3fc11112, 0xe0a6b45f } else .{ 0xe0a6b45f, 0x3fc11112 }; //    .
// constants
pub const mfftnhf: [2]u32 = if (endian == .big) .{ 0xc02f0000, 0x00000000 } else .{ 0x00000000, 0xc02f0000 }; //-15.5
pub const g1: [2]u32 = if (endian == .big) .{ 0x3e4b096c, 0x00000000 } else .{ 0x00000000, 0x3e4b096c }; // 1.259e-8
pub const g2: [2]u32 = if (endian == .big) .{ 0x3faf212d, 0x00000000 } else .{ 0x00000000, 0x3faf212d }; // 0.0608
pub const g3: [2]u32 = if (endian == .big) .{ 0x3fe92f1a, 0x00000000 } else .{ 0x00000000, 0x3fe92f1a }; // 0.787
pub const g4: [2]u32 = if (endian == .big) .{ 0x40390000, 0x00000000 } else .{ 0x00000000, 0x40390000 }; // 25.0
pub const g5: [2]u32 = if (endian == .big) .{ 0x4197d784, 0x00000000 } else .{ 0x00000000, 0x4197d784 }; // 1e8
pub const gy2: [2]u32 = if (endian == .big) .{ 0x3faf212d, 0x00000000 } else .{ 0x00000000, 0x3faf212d }; // 0.0608
pub const mp1: [2]u32 = if (endian == .big) .{ 0x3ff921fb, 0x58000000 } else .{ 0x58000000, 0x3ff921fb };
pub const mp2: [2]u32 = if (endian == .big) .{ 0xbe4dde97, 0x3c000000 } else .{ 0x3c000000, 0xbe4dde97 };
pub const mp3: [2]u32 = if (endian == .big) .{ 0xbc8cb3b3, 0x99d747f2 } else .{ 0x99d747f2, 0xbc8cb3b3 };
pub const pp3: [2]u32 = if (endian == .big) .{ 0xbc8cb3b3, 0x98000000 } else .{ 0x98000000, 0xbc8cb3b3 };
pub const pp4: [2]u32 = if (endian == .big) .{ 0xbacd747f, 0x23e32ed7 } else .{ 0x23e32ed7, 0xbacd747f };
pub const hpinv: [2]u32 = if (endian == .big) .{ 0x3fe45f30, 0x6dc9c883 } else .{ 0x6dc9c883, 0x3fe45f30 };
pub const toint: [2]u32 = if (endian == .big) .{ 0x43380000, 0x00000000 } else .{ 0x00000000, 0x43380000 };

pub const xfg = @import("utan_tbl.zig").xfg;
