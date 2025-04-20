// tab[i] = uint(2^(i/32)) - (i << 52-BITS)
// used for computing 2^(k/32) for an int |k| < 150 * 32 as
// double(tab[k%32] + (k << 52-BITS))
pub const T_32: [32]u64 = .{
    0x3ff0000000000000, 0x3fefd9b0d3158574, 0x3fefb5586cf9890f, 0x3fef9301d0125b51,
    0x3fef72b83c7d517b, 0x3fef54873168b9aa, 0x3fef387a6e756238, 0x3fef1e9df51fdee1,
    0x3fef06fe0a31b715, 0x3feef1a7373aa9cb, 0x3feedea64c123422, 0x3feece086061892d,
    0x3feebfdad5362a27, 0x3feeb42b569d4f82, 0x3feeab07dd485429, 0x3feea47eb03a5585,
    0x3feea09e667f3bcd, 0x3fee9f75e8ec5f74, 0x3feea11473eb0187, 0x3feea589994cce13,
    0x3feeace5422aa0db, 0x3feeb737b0cdc5e5, 0x3feec49182a3f090, 0x3feed503b23e255d,
    0x3feee89f995ad3ad, 0x3feeff76f2fb5e47, 0x3fef199bdd85529c, 0x3fef3720dcef9069,
    0x3fef5818dcfba487, 0x3fef7c97337b9b5f, 0x3fefa4afa2a490da, 0x3fefd0765b6e4540,
};
pub const shift_scaled_32: f64 = @as(f64, 0x1.8p+52) / 32.0;
pub const poly: [3]f64 = .{ 0x1.c6af84b912394p-5, 0x1.ebfce50fac4f3p-3, 0x1.62e42ff0c52d6p-1 };
pub const SHIFT_32: f64 = 0x1.8p+52;
pub const InvLn2N_32: f64 = 0x1.71547652b82fep+0 * 32.0;
pub const C_32: [3]f64 = .{ 0x1.c6af84b912394p-5 / 32.0 / 32.0 / 32.0, 0x1.ebfce50fac4f3p-3 / 32.0 / 32.0, 0x1.62e42ff0c52d6p-1 / 32.0 };
