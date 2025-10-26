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
    var size = 0;
    if (X == comptime_int) {
        if (x == 0)
            size = 1
        else
            size = std.math.log2(int.abs(x)) / 32 + 1;
    } else {
        switch (@typeInfo(X).int.bits) {
            8, 16, 32 => size = 1,
            64 => size = 2,
            128 => size = 4,
            else => unreachable,
        }
    }

    break :t struct { integer.Integer, [size]u32 };
} {
    const X: type = @TypeOf(x);

    comptime if (types.numericType(X) != .int)
        @compileError("int.asInteger requires x to be an int type, got " ++ @typeName(X));

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
            limbs[i] = @truncate(uvalue & 0xFFFFFFFF);
            uvalue >>= 32;
        }

        //std.debug.print("limbs = {any}\n", .{limbs[0..size]});

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
    const size: u32 = if (bits <= 32) 1 else bits / 32;

    var limbs: [size]u32 = undefined;

    const uvalue: std.meta.Int(.unsigned, bits) = types.scast(std.meta.Int(.unsigned, bits), int.abs(x));

    var actual_size: u32 = 0;
    switch (bits) {
        8, 16, 32 => {
            if (uvalue != 0) {
                limbs[0] = uvalue;
                actual_size = 1;
            }
        },
        64, 128 => {
            if (x != 0) {
                var chunks: [bits / 32]u32 = @bitCast(uvalue);
                if (comptime builtin.cpu.arch.endian() == .big) {
                    // Reverse
                    var i: u32 = 0;
                    while (i < (bits / 32) / 2) : (i += 1) {
                        const temp: u32 = chunks[i];
                        chunks[i] = chunks[bits / 32 - 1 - i];
                        chunks[bits / 32 - 1 - i] = temp;
                    }
                }

                var i: u32 = 0;
                while (i < bits / 32 and chunks[i] != 0) : (i += 1) {
                    limbs[i] = chunks[i];
                    actual_size += 1;
                }
            }
        },
        else => unreachable,
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
