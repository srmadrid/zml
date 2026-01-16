const std = @import("std");
const builtin = @import("builtin");

const types = @import("../types.zig");
const int = @import("../int.zig");
const integer = @import("../integer.zig");

/// Only for internal use.
///
/// Converts an int type to a `Integer` as a view. After calling, you must execute:
/// ```zig
/// var i = asInteger(x);
/// i[0].limbs = &i[1];
/// ```
pub fn asInteger(x: anytype) t: {
    const X: type = @TypeOf(x);

    if (!types.isNumeric(X) or types.numericType(X) != .int)
        @compileError("zml.int.asInteger: x must be an int, got \n\tx: " ++ @typeName(X) ++ "\n");

    var size = 0;
    if (X == comptime_int) {
        if (x == 0)
            size = 1
        else
            size = std.math.log2(int.abs(x)) / 32 + 1;
    } else {
        size = int.max(1, @typeInfo(X).int.bits / 32);
    }

    break :t struct { integer.Integer, [size]u32 };
} {
    const X: type = @TypeOf(x);

    if (comptime X == comptime_int) {
        if (x == 0) {
            return .{
                .{
                    .limbs = undefined,
                    .size = 1,
                    ._llen = 1,
                    .positive = true,
                    .flags = .{ .owns_data = false, .writable = false },
                },
                .{0},
            };
        }

        comptime var uvalue: comptime_int = int.abs(x);
        const size: u32 = types.scast(u32, std.math.log2(uvalue) / 32 + 1);

        var limbs: [size]u32 = undefined;

        comptime var i: u32 = 0;
        inline while (uvalue != 0) : (i += 1) {
            limbs[i] = @truncate(uvalue & 0xffffffff);
            uvalue >>= 32;
        }

        return .{
            .{
                .limbs = undefined,
                .size = size,
                ._llen = size,
                .positive = x >= 0,
                .flags = .{ .owns_data = false, .writable = false },
            },
            limbs,
        };
    }

    const bits: u16 = @typeInfo(X).int.bits;
    const size: u32 = int.max(1, bits / 32);

    var limbs: [size]u32 = undefined;

    var uvalue: std.meta.Int(.unsigned, bits) = types.scast(std.meta.Int(.unsigned, bits), int.abs(x));

    var actual_size: u32 = 0;
    if (bits <= 32) {
        if (uvalue != 0) {
            limbs[0] = uvalue;
            actual_size = 1;
        }
    } else {
        if (x != 0) {
            var chunks: [size]u32 = undefined;
            while (uvalue != 0) {
                chunks[actual_size] = @truncate(uvalue & 0xffffffff);
                uvalue >>= 32;
                actual_size += 1;
            }

            if (comptime builtin.cpu.arch.endian() == .big) {
                // Reverse
                var i: u32 = 0;
                while (i < size / 2) : (i += 1) {
                    const temp: u32 = chunks[i];
                    chunks[i] = chunks[size - 1 - i];
                    chunks[size - 1 - i] = temp;
                }
            }

            var i: u32 = 0;
            while (i < size and chunks[i] != 0) : (i += 1) {
                limbs[i] = chunks[i];
                actual_size += 1;
            }
        }
    }

    return .{
        .{
            .limbs = undefined,
            .size = actual_size,
            ._llen = size,
            .positive = x >= 0,
            .flags = .{ .owns_data = false, .writable = false },
        },
        limbs,
    };
}
