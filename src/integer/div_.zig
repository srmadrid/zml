const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Performs in-place division between two operands of any numeric type in
/// `Integer` precision. Float, rational or real types are truncated towards
/// zero, and for cfloat or complex types, only the real part is considered.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// Signature
/// ---------
/// ```zig
/// fn div_(allocator: std.mem.Allocator, o: *Integer, x: X, y: Y, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `allocator` (`std.mem.Allocator`):
/// The allocator to use for memory allocations. Must be the same allocator used
/// to initialize `o`.
///
/// `o` (`*Integer`):
/// A pointer to the output operand where the result will be stored.
///
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// Returns
/// -------
/// `void`
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails.
///
/// `integer.Error.NotWritable`:
/// If the output operand `o` is not writable.
///
/// `integer.Error.ZeroDivision`:
/// If `y` is zero.
pub fn div_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("integer.div_ requires x and y to be numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return integer.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .complex => switch (comptime types.numericType(Y)) {
            .complex => return div_(allocator, o, x.re, y.re),
            .real => @compileError("integer.div_ not implemented for Complex + Real yet"),
            .rational => return div_(allocator, o, x.re, y),
            .integer => return div_(allocator, o, x.re, y),
            .cfloat => return div_(allocator, o, x.re, y.re),
            .float => return div_(allocator, o, x.re, y),
            .int => return div_(allocator, o, x.re, y),
            .bool => return div_(allocator, o, x.re, y),
        },
        .real => @compileError("integer.div_ not implemented for Real yet"),
        .rational => switch (comptime types.numericType(Y)) {
            .complex => return div_(allocator, o, x, y.re),
            .real => @compileError("integer.div_ not implemented for Rational + Real yet"),
            .rational => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return div_(allocator, o, tx, ty);
            },
            .integer => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                return div_(allocator, o, tx, y);
            },
            .cfloat => return div_(allocator, o, x, y.re),
            .float => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return div_(allocator, o, tx, ty[0]);
            },
            .int => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return div_(allocator, o, tx, ty[0]);
            },
            .bool => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                return div_(allocator, o, tx, types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .integer => switch (comptime types.numericType(Y)) {
            .complex => return div_(allocator, o, x, y.re),
            .real => @compileError("integer.div_ not implemented for Integer + Real yet"),
            .rational => {
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return div_(allocator, o, x, ty);
            },
            .integer => {
                if (y.size == 0)
                    return integer.Error.ZeroDivision;

                if (integer.lt(integer.abs(null, x) catch unreachable, integer.abs(null, y) catch unreachable))
                    return o.set(allocator, 0);

                if (y.size == 1) {
                    try o.reserve(allocator, x.size);
                    var carry: u64 = 0;
                    var i: u32 = x.size;
                    while (i > 0) : (i -= 1) {
                        const limb: u64 = (carry << 32) | types.scast(u64, x.limbs[i - 1]);
                        const q: u32 = @truncate(limb / types.scast(u64, y.limbs[0]));
                        carry = limb % types.scast(u64, y.limbs[0]);
                        o.limbs[i - 1] = q;
                    }

                    o.size = x.size;
                    o.truncate();
                    o.positive = x.positive == y.positive;
                    return;
                }

                // x.size and y.size > 1
                const n: u32 = y.size;
                const m: u32 = x.size - n;

                // Ensure o has enough capacity for the left shift
                var r: Integer = try .init(allocator, x.size + 1);
                try r.set(allocator, x);
                defer r.deinit(allocator);

                const s: u32 = @clz(y.limbs[n - 1]);

                // We left shift so y.limbs[n - 1] & (1 << 31) != 0
                // So, no need to allocate an extra limb for v
                var v: Integer = try y.copy(allocator);
                defer v.deinit(allocator);

                if (s > 0) {
                    _ = shiftLeftInPlace(v.limbs[0..v.size], null, @intCast(s));
                    if (shiftLeftInPlace(r.limbs[0..r.size], &r.limbs[r.size], @intCast(s)))
                        r.size += 1;
                }

                try o.reserve(allocator, m + 1);

                var j: u32 = m;
                while (true) : (j -= 1) {
                    const qhat: u32 = estimate_qhat(
                        if (j + n < r.size) r.limbs[j + n] else 0,
                        if (j + n - 1 < r.size) r.limbs[j + n - 1] else 0,
                        if (j + n - 2 < r.size) r.limbs[j + n - 2] else 0,
                        v.limbs[n - 1],
                        v.limbs[n - 2],
                    );

                    o.limbs[j] = qhat;

                    if (mul_sub_(r.limbs[j .. j + n + 1], v.limbs[0..n], qhat)) {
                        o.limbs[j] -= 1;
                        add_back(r.limbs[j .. j + n + 1], v.limbs[0..n]);
                    }

                    if (j == 0) break;
                }

                o.size = m + 1;
                o.truncate();

                if (s > 0) shiftRightInPlace(r.limbs[0..n], @intCast(s));
                if (r.limbs[r.size - 1] == 0 and r.size > 1) r.size -= 1;

                o.positive = x.positive == y.positive;
            },
            .cfloat => return div_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return div_(allocator, o, x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return div_(allocator, o, x, ty[0]);
            },
            .bool => {
                return div_(allocator, o, x, types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .complex => return div_(allocator, o, x.re, y.re),
            .real => @compileError("integer.div_ not implemented for CFloat + Real yet"),
            .rational => return div_(allocator, o, x.re, y),
            .integer => return div_(allocator, o, x.re, y),
            .cfloat => return div_(allocator, o, x.re, y.re),
            .float => return div_(allocator, o, x.re, y),
            .int => return div_(allocator, o, x.re, y),
            .bool => return div_(allocator, o, x.re, y),
        },
        .float => switch (comptime types.numericType(Y)) {
            .complex => return div_(allocator, o, x, y.re),
            .real => @compileError("integer.div_ not implemented for Float + Real yet"),
            .rational => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return div_(allocator, o, tx[0], ty);
            },
            .integer => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return div_(allocator, o, tx[0], y);
            },
            .cfloat => return div_(allocator, o, x, y.re),
            .float => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return div_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return div_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return div_(allocator, o, tx[0], types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .complex => return div_(allocator, o, x, y.re),
            .real => @compileError("integer.div_ not implemented for Int + Real yet"),
            .rational => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return div_(allocator, o, tx[0], ty);
            },
            .integer => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return div_(allocator, o, tx[0], y);
            },
            .cfloat => return div_(allocator, o, x, y.re),
            .float => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return div_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return div_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return div_(allocator, o, tx[0], types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .complex => return div_(allocator, o, x, y.re),
            .real => @compileError("integer.div_ not implemented for Bool + Real yet"),
            .rational => {
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return div_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, ty);
            },
            .integer => {
                return div_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, y);
            },
            .cfloat => return div_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return div_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, ty[0]);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return div_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, ty[0]);
            },
            .bool => {
                return div_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    types.cast(Integer, y, .{}) catch unreachable,
                );
            },
        },
    }
}

pub fn shiftLeftInPlace(limbs: []u32, next: ?*u32, s: u5) bool {
    if (s == 0 or limbs.len == 0) return false;

    var carry: u32 = 0;
    var i: u32 = 0;
    while (i < limbs.len) : (i += 1) {
        const new_carry: u32 = limbs[i] >> @as(u5, @intCast(32 - @as(u6, @intCast(s))));
        limbs[i] = (limbs[i] << s) | carry;
        carry = new_carry;
    }

    if (carry != 0) {
        // Extend by one limb if necessary
        if (next) |n| {
            n.* = carry;
        }

        return true;
    }

    return false;
}

pub fn shiftRightInPlace(limbs: []u32, s: u5) void {
    if (s == 0 or limbs.len == 0) return;

    var carry: u32 = 0;
    var i: u32 = types.scast(u32, limbs.len);
    while (i > 0) : (i -= 1) {
        const limb: u32 = limbs[i - 1];
        const new_carry: u32 = limb << @as(u5, @intCast(32 - @as(u6, @intCast(s))));
        limbs[i - 1] = (limb >> s) | carry;
        carry = new_carry;
    }
}

fn estimate_qhat(
    u_high: u32, // u[j+n]
    u_next: u32, // u[j+n−1]
    u_next2: u32, // u[j+n−2]
    v_high: u32, // v[n−1]
    v_next: u32, // v[n−2]
) u32 {
    const b: u64 = 1 << 32;
    const num: u64 = (types.scast(u64, u_high) << 32) | u_next;
    var qhat: u64 = num / types.scast(u64, v_high);
    var rhat: u64 = num % types.scast(u64, v_high);

    var i: usize = 0;
    while (i < 2) : (i += 1) {
        if (qhat >= b or (qhat * types.scast(u64, v_next) > (rhat << 32) + types.scast(u64, u_next2))) {
            qhat -= 1;
            rhat += types.scast(u64, v_high);
            if (rhat >= b) break;
        } else break;
    }

    return @truncate(qhat);
}

fn mul_sub_(u_part: []u32, v: []const u32, qhat: u32) bool {
    var carry: u64 = 0;
    var borrow: u64 = 0;
    const n = v.len;

    var i: usize = 0;
    while (i < n) : (i += 1) {
        const prod: u64 = types.scast(u64, v[i]) * qhat + carry;
        carry = prod >> 32;
        const prod_low: u64 = prod & 0xFFFFFFFF;

        // Compute u_part[i] - prod_low - borrow
        const u_val: u64 = @as(u64, u_part[i]);
        const sub_val: u64 = prod_low + borrow;

        if (u_val >= sub_val) {
            u_part[i] = @truncate(u_val - sub_val);
            borrow = 0;
        } else {
            u_part[i] = @truncate(u_val + (1 << 32) - sub_val);
            borrow = 1;
        }
    }

    const high_res = @subWithOverflow(@as(u64, u_part[n]), carry + borrow);
    u_part[n] = @truncate(high_res[0]);

    return high_res[1] != 0;
}

fn add_back(u_part: []u32, v: []const u32) void {
    var carry: u64 = 0;
    const n = v.len;
    var i: usize = 0;

    while (i < n) : (i += 1) {
        const sum: u64 = @as(u64, u_part[i]) + @as(u64, v[i]) + carry;
        u_part[i] = @truncate(sum);
        carry = sum >> 32;
    }

    var k = n;
    while (carry != 0 and k < u_part.len) : (k += 1) {
        const sum: u64 = @as(u64, u_part[k]) + carry;
        u_part[k] = @truncate(sum);
        carry = sum >> 32;
    }
}
