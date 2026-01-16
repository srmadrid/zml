const std = @import("std");
const builtin = @import("builtin");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

/// Only for internal use.
///
/// Converts a float type to a `Rational` as a view. After calling, you must execute:
/// ```zig
/// var r = asRational(x);
/// r[0].num.limbs = &r[1][0];
/// r[0].den.limbs = &r[1][1];
/// ```
pub fn asRational(x: anytype) !t: {
    const X: type = @TypeOf(x);

    if (!types.isNumeric(X) or types.numericType(X) != .float)
        @compileError("zml.float.asRational: x must be a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    var size = 0;
    if (X == comptime_float) {
        size = 512;
    } else {
        switch (@typeInfo(X).float.bits) {
            16 => size = 1,
            32 => size = 4,
            64 => size = 32,
            80 => size = 512,
            128 => size = 512,
            else => unreachable,
        }
    }

    break :t struct { rational.Rational, [2][size]u32 };
} {
    const X: type = @TypeOf(x);

    const v: if (X == comptime_float) f128 else X = x;

    if (!std.math.isFinite(v))
        return float.Error.NotFinite;

    const bits: u16 = @typeInfo(@TypeOf(v)).float.bits;
    const size: u32 = switch (bits) {
        16 => 1,
        32 => 4,
        64 => 32,
        80 => 512,
        128 => 512,
        else => unreachable,
    };

    var nlimbs: [size]u32 = undefined;
    var dlimbs: [size]u32 = undefined;

    if (v == 0) {
        dlimbs[0] = 1;

        return .{
            .{
                .num = .{
                    .limbs = undefined,
                    .size = 0,
                    ._llen = 1,
                    .positive = true,
                    .flags = .{ .owns_data = false, .writable = false },
                },
                .den = .{
                    .limbs = undefined,
                    .size = 1,
                    ._llen = 1,
                    .positive = true,
                    .flags = .{ .owns_data = false, .writable = false },
                },
                .flags = .{ .owns_data = false, .writable = false },
            },
            .{ nlimbs, dlimbs },
        };
    }

    const U: type = switch (bits) {
        16 => u16,
        32 => u32,
        64 => u64,
        80 => u80,
        128 => u128,
        else => unreachable,
    };

    const E: type = switch (bits) {
        16 => i6,
        32 => i9,
        64 => i12,
        80 => i16,
        128 => i16,
        else => unreachable,
    };
    const einfo = @typeInfo(E);
    const exp_bias: E = (1 << (einfo.int.bits - 2)) - 1;
    const exp_mask: U = (1 << einfo.int.bits - 1) - 1;

    const F: type = switch (bits) {
        16 => u10,
        32 => u23,
        64 => u52,
        80 => u64,
        128 => u112,
        else => unreachable,
    };
    const finfo = @typeInfo(F);
    const frac_mask: U = (1 << finfo.int.bits) - 1;

    const uvalue: U = @bitCast(v);
    const sign_bit: bool = (uvalue >> (bits - 1)) != 0;
    const raw_exp: U = (uvalue >> finfo.int.bits) & exp_mask;
    var mantissa: U = uvalue & frac_mask;
    mantissa |= (@as(U, 1) << finfo.int.bits); // implicit leading one

    var exponent: E = @as(E, @intCast(raw_exp)) - exp_bias;

    if (raw_exp == 0) {
        if (mantissa == 0)
            exponent = 0
        else
            exponent = -exp_bias + 1;
    } else {
        mantissa |= (@as(U, 1) << finfo.int.bits);
    }

    exponent -= @as(E, @intCast(finfo.int.bits));

    var nsize: u32 = 0;
    var mtmp: if (U == u16 or U == u32) u64 else U = @intCast(mantissa);
    while (mtmp != 0) : (mtmp >>= 32) {
        nlimbs[nsize] = @truncate(mtmp & 0xFFFFFFFF);
        nsize += 1;
    }

    var dsize: u32 = 1;
    dlimbs[0] = 1;

    if (exponent > 0) {
        var shift: u32 = @intCast(exponent);

        // Shift by whole limbs (32-bit steps) first
        while (shift >= 32) : (shift -= 32) {
            // shift entire limb array left by 32 bits (insert zero limb at start)
            var i: u32 = nsize;
            while (i > 0) : (i -= 1) {
                nlimbs[i] = nlimbs[i - 1];
            }
            nlimbs[0] = 0;
            nsize += 1;
        }

        // Remaining partial bits (0..31)
        if (shift > 0) {
            var extra: u32 = 0;
            const carry = @import("../integer/div_.zig").shiftLeftInPlace(nlimbs[0..nsize], &extra, @intCast(shift & 31));
            if (carry) {
                nlimbs[nsize] = extra;
                nsize += 1;
            }
        }
    } else if (exponent < 0) {
        var shift: u32 = @intCast(-exponent);

        // Shift denominator by powers of two
        while (shift >= 32) : (shift -= 32) {
            var i: u32 = dsize;
            while (i > 0) : (i -= 1) {
                dlimbs[i] = dlimbs[i - 1];
            }
            dlimbs[0] = 0;
            dsize += 1;
        }

        if (shift > 0) {
            var extra: u32 = 0;
            const carry = @import("../integer//div_.zig").shiftLeftInPlace(dlimbs[0..dsize], &extra, @intCast(shift & 31));
            if (carry) {
                dlimbs[dsize] = extra;
                dsize += 1;
            }
        }
    }

    return .{
        .{
            .num = .{
                .limbs = undefined,
                .size = nsize,
                ._llen = size,
                .positive = !sign_bit,
                .flags = .{ .owns_data = false, .writable = false },
            },
            .den = .{
                .limbs = undefined,
                .size = dsize,
                ._llen = size,
                .positive = true,
                .flags = .{ .owns_data = false, .writable = false },
            },
            .flags = .{ .owns_data = false, .writable = false },
        },
        .{ nlimbs, dlimbs },
    };
}
