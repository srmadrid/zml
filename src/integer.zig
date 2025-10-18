const std = @import("std");
const builtin = @import("builtin");

const types = @import("types.zig");
const constants = @import("constants.zig");
const int = @import("int.zig");
const float = @import("float.zig");
const rational = @import("rational.zig");
const Rational = rational.Rational;

pub const Integer = struct {
    limbs: [*]u32,
    size: u32,
    _llen: u32,
    positive: bool,
    flags: Flags,

    pub const empty: Integer = .{
        .limbs = &.{},
        .size = 0,
        ._llen = 0,
        .positive = true,
        .flags = .{ .owns_data = false, .writable = false },
    };

    pub fn init(allocator: std.mem.Allocator, size: u32) !Integer {
        // if (size == 0)
        //     return Error.ZeroSize;

        return .{
            .limbs = (try allocator.alloc(u32, size)).ptr,
            .size = 0,
            ._llen = size,
            .positive = true,
            .flags = .{ .owns_data = true, .writable = true },
        };
    }

    pub fn initSet(allocator: std.mem.Allocator, value: anytype) !Integer {
        var integer = try Integer.init(allocator, 0);

        try integer.set(allocator, value);
        return integer;
    }

    pub fn deinit(self: *Integer, allocator: ?std.mem.Allocator) void {
        if (self.flags.owns_data) {
            allocator.?.free(self.limbs[0..self._llen]);
        }

        self.* = undefined;
    }

    pub fn reserve(self: *Integer, allocator: std.mem.Allocator, new_size: u32) !void {
        if (self.flags.owns_data == false)
            return Error.DataNotOwned;

        if (new_size > self._llen) {
            self.limbs = (try allocator.realloc(self.limbs[0..self._llen], new_size)).ptr;
            self._llen = new_size;
        }
    }

    pub fn trim(self: *Integer, allocator: std.mem.Allocator) !void {
        if (self.flags.owns_data == false or self.size == self._llen)
            return;

        self.limbs = (try allocator.realloc(self.limbs[0..self._llen], self.size)).ptr;
        self._llen = self.size;
    }

    pub fn trimSize(self: *Integer) void {
        while (self.size > 0 and self.limbs[self.size - 1] == 0) {
            self.size -= 1;
        }

        if (self.size == 0) self.positive = true;
    }

    pub fn set(self: *Integer, allocator: std.mem.Allocator, value: anytype) !void {
        const V: type = @TypeOf(value);

        if (!self.flags.writable)
            return Error.NotWritable;

        if (comptime types.isNumeric(V)) {
            switch (comptime types.numericType(V)) {
                .bool => {
                    try self.reserve(allocator, 1);

                    if (self) {
                        self.limbs[0] = 1;
                        self.size = 1;
                        self.positive = true;
                    } else {
                        self.size = 0;
                        self.positive = true;
                    }
                },
                .int => {
                    if (value == 0) {
                        try self.reserve(allocator, 1);

                        self.size = 0;
                        self.positive = true;
                        return;
                    }

                    if (comptime V == comptime_int) {
                        comptime var abs_value: comptime_int = int.abs(value);
                        const size: u32 = types.scast(u32, std.math.log2(abs_value) / 32 + 1);

                        if (self.flags.owns_data) {
                            try self.reserve(allocator, size);
                        } else if (self._llen < size) {
                            return Error.DataNotOwned;
                        }

                        comptime var i: u32 = 0;
                        inline while (abs_value != 0) : (i += 1) {
                            self.limbs[i] = @truncate(abs_value & 0xFFFFFFFF);
                            self.size += 1;
                            abs_value >>= 32;
                        }

                        if (value < 0) self.positive = false;

                        return;
                    }

                    const bits: u16 = @typeInfo(V).int.bits;
                    const size: u32 = if (bits <= 32) 1 else bits / 32;

                    if (self.flags.owns_data) {
                        try self.reserve(allocator, size);
                    } else if (self._llen < size) {
                        return Error.DataNotOwned;
                    }

                    if (value < 0) self.positive = false;

                    switch (bits) {
                        8 => {
                            const uvalue: u8 = if (comptime @typeInfo(V).int.signedness == .signed)
                                value & 0x7F
                            else
                                value;

                            if (uvalue != 0) {
                                self.limbs[0] = uvalue;
                                self.size = 1;
                            }
                        },
                        16 => {
                            const uvalue: u16 = if (comptime @typeInfo(V).int.signedness == .signed)
                                value & 0x7FFF
                            else
                                value;

                            if (uvalue != 0) {
                                self.limbs[0] = uvalue;
                                self.size = 1;
                            }
                        },
                        32, 64, 128 => {
                            if (value != 0) {
                                var chunks: [bits / 32]u32 = @bitCast(value);
                                if (comptime builtin.cpu.arch.endian() == .big) {
                                    // Reverse
                                    var i: u32 = 0;
                                    while (i < (bits / 32) / 2) : (i += 1) {
                                        const temp: u32 = chunks[i];
                                        chunks[i] = chunks[bits / 32 - 1 - i];
                                        chunks[bits / 32 - 1 - i] = temp;
                                    }
                                }

                                if (comptime @typeInfo(V).int.signedness == .signed)
                                    chunks[bits / 32 - 1] &= 0x7FFFFFFF; // Remove sign

                                var i: u32 = 0;
                                while (i < bits / 32) : (i += 1) {
                                    self.limbs[i] = chunks[i];
                                }

                                self.size = bits / 32;
                            }
                        },
                        else => unreachable,
                    }
                },
                .float => {
                    const v: V = float.trunc(value);

                    if (v == 0) {
                        try self.reserve(allocator, 1);

                        self.size = 0;
                        self.positive = true;
                        return;
                    }

                    if (comptime V == comptime_float) {
                        return self.set(allocator, comptime types.scast(comptime_int, v));
                    }

                    if (!std.math.isFinite(v))
                        return Error.NonInteger;

                    const bits: u16 = @typeInfo(V).float.bits;

                    switch (bits) {
                        16 => {
                            const size: u32 = 1; // Max size for f16 is 1 limb (u32)

                            if (self.flags.owns_data) {
                                try self.reserve(allocator, size);
                            } else if (self._llen < size) {
                                return Error.DataNotOwned;
                            }

                            if (value < 0) self.positive = false;

                            const uvalue: u16 = @bitCast(v);
                            const exponent: i16 = types.scast(i16, (uvalue & 0x7C00) >> 10) - 15;
                            var mantissa: u32 = types.scast(u32, uvalue & 0x3FF);

                            if (exponent == -15) {
                                return; // Zero or subnormal
                            } else {
                                mantissa |= 0x400;
                            }

                            if (exponent > 10) {
                                mantissa <<= @as(u5, @intCast(exponent - 10));
                            } else {
                                mantissa >>= @as(u5, @intCast(10 - exponent));
                            }

                            // Mantissa is now the integer value
                            if (mantissa != 0) {
                                self.limbs[0] = mantissa;
                                self.size = 1;
                            }
                        },
                        32 => {
                            const uvalue: u32 = @bitCast(v);
                            const exponent: i32 = types.scast(i32, (uvalue & 0x7F800000) >> 23) - 127;
                            var mantissa: u160 = @intCast(uvalue & 0x7FFFFF);

                            const size: u32 = if (exponent <= 8) 1 else (types.scast(u32, exponent) / 32 + 1);
                            if (self.flags.owns_data) {
                                try self.reserve(allocator, size);
                            } else if (self._llen < size) {
                                return Error.DataNotOwned;
                            }

                            if (value < 0) self.positive = false;

                            if (exponent == -127) {
                                return; // Zero or subnormal
                            } else {
                                mantissa |= 0x800000;
                            }

                            if (exponent > 23) {
                                mantissa <<= @as(u7, @intCast(exponent - 23));
                            } else if (exponent >= 0) {
                                mantissa >>= @as(u7, @intCast(23 - exponent));
                            } else {
                                mantissa = 0;
                            }

                            // Mantissa is now the integer value
                            while (mantissa != 0) {
                                self.limbs[self.size] = @truncate(mantissa & 0xFFFFFFFF);
                                self.size += 1;
                                mantissa >>= 32;
                            }
                        },
                        64 => {
                            const uvalue: u64 = @bitCast(v);
                            const exponent: i32 = types.scast(i32, (uvalue & 0x7FF0000000000000) >> 52) - 1023;
                            var mantissa: u1088 = @intCast(uvalue & 0xFFFFFFFFFFFFF);

                            const size: u32 = if (exponent <= -21) 1 else (types.scast(u32, exponent) / 32 + 1);
                            if (self.flags.owns_data) {
                                try self.reserve(allocator, size);
                            } else if (self._llen < size) {
                                return Error.DataNotOwned;
                            }

                            if (value < 0) self.positive = false;

                            if (exponent == -1023) {
                                return; // Zero or subnormal
                            } else {
                                mantissa |= 0x10000000000000;
                            }

                            if (exponent > 52) {
                                mantissa <<= @as(u11, @intCast(exponent - 52));
                            } else if (exponent >= 0) {
                                mantissa >>= @as(u11, @intCast(52 - exponent));
                            } else {
                                mantissa = 0;
                            }

                            // Mantissa is now the integer value
                            while (mantissa != 0) {
                                self.limbs[self.size] = @truncate(mantissa & 0xFFFFFFFF);
                                self.size += 1;
                                mantissa >>= 32;
                            }
                        },
                        80 => {
                            return self.set(allocator, types.scast(f128, value));
                        },
                        128 => {
                            const uvalue: u128 = @bitCast(v);
                            const exponent: i32 = types.scast(i32, (uvalue & 0x7FFF0000000000000000000000000000) >> 112) - 16383;
                            var mantissa: u16512 = @intCast(uvalue & 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFF);

                            const size: u32 = if (exponent <= -81) 1 else (types.scast(u32, exponent) / 32 + 1);
                            if (self.flags.owns_data) {
                                try self.reserve(allocator, size);
                            } else if (self._llen < size) {
                                return Error.DataNotOwned;
                            }

                            if (value < 0) self.positive = false;

                            if (exponent == -16383) {
                                return; // Zero or subnormal
                            } else {
                                mantissa |= 0x10000000000000000000000000000;
                            }

                            if (exponent > 112) {
                                mantissa <<= @as(u15, @intCast(exponent - 112));
                            } else if (exponent >= 0) {
                                mantissa >>= @as(u15, @intCast(112 - exponent));
                            } else {
                                mantissa = 0;
                            }

                            // Mantissa is now the integer value
                            while (mantissa != 0) {
                                self.limbs[self.size] = @truncate(mantissa & 0xFFFFFFFF);
                                self.size += 1;
                                mantissa >>= 32;
                            }
                        },
                        else => unreachable,
                    }
                },
                .cfloat => {
                    return self.set(allocator, value.re);
                },
                .integer => {
                    if (self.flags.owns_data) {
                        try self.reserve(allocator, value.size);
                    } else if (self._llen < value.size) {
                        return Error.DataNotOwned;
                    }

                    var i: u32 = 0;
                    while (i < value.size) : (i += 1) {
                        self.limbs[i] = value.limbs[i];
                    }
                    self.size = value.size;
                    self.positive = value.positive;

                    return;
                },
                .rational => {
                    if (eq(value.den, 1, .{}) catch unreachable) {
                        return self.set(allocator, value.num);
                    } else {
                        return div_(allocator, self, value.num, value.den);
                    }
                },
                .real => {},
                .complex => {
                    return self.set(allocator, value.re);
                },
                .expression => {},
            }
        } else if (comptime V == []const u8 or V == []u8) {} else @compileError("Value must be a numeric type or a string");
    }

    pub fn copy(self: *const Integer, allocator: std.mem.Allocator) !Integer {
        var result: Integer = try .init(allocator, int.max(1, self.size));

        var i: u32 = 0;
        while (i < self.size) : (i += 1) {
            result.limbs[i] = self.limbs[i];
        }

        result.size = self.size;
        result.positive = self.positive;

        return result;
    }

    pub fn asRational(self: *const Integer) Rational {
        var num: Integer = self.*;
        num.flags.owns_data = false;

        return .{
            .num = num,
            .den = constants.one(Integer, .{}) catch unreachable,
            .flags = .{ .owns_data = false, .writable = false },
        };
    }
};

// Arithmetic operations
pub const add = @import("integer/add.zig").add;
pub const add_ = @import("integer/add_.zig").add_;
pub const sub = @import("integer/sub.zig").sub;
pub const sub_ = @import("integer/sub_.zig").sub_;
pub const mul = @import("integer/mul.zig").mul;
pub const mul_ = @import("integer/mul_.zig").mul_;
pub const div = @import("integer/div.zig").div;
pub const div_ = @import("integer/div_.zig").div_;

// Comparison operations
pub const cmp = @import("integer/cmp.zig").cmp;
pub const eq = @import("integer/eq.zig").eq;
pub const ne = @import("integer/ne.zig").ne;
pub const lt = @import("integer/lt.zig").lt;
pub const le = @import("integer/le.zig").le;
pub const gt = @import("integer/gt.zig").gt;
pub const ge = @import("integer/ge.zig").ge;

// Basic operations
pub const abs = @import("integer/abs.zig").abs;
pub const neg = @import("integer/neg.zig").neg;

pub const gcd = @import("integer/gcd.zig").gcd;

pub const Error = error{
    ZeroSize,
    ZeroDivision,
    NonInteger,
    NotWritable,
    DataNotOwned,
};

pub const Flags = packed struct {
    owns_data: bool = true,
    writable: bool = true,
};
