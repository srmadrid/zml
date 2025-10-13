const std = @import("std");

const types = @import("../types.zig");
const constants = @import("../constants.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn div_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("integer.add_ requires x and y to be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                if (x.size == 0) return integer.Error.ZeroDivision;

                if (integer.lt(integer.abs(null, x) catch unreachable, integer.abs(null, y) catch unreachable))
                    return o.set(allocator, 0);

                const n: u32 = y.size;
                const m: u32 = x.size - n;

                try o.reserve(allocator, m + 1);
                var r: Integer = try x.copy(allocator);
                defer r.deinit(allocator);

                const s: u32 = @clz(y.limbs[n - 1]);
                var v: Integer = try y.copy(allocator);
                defer v.deinit(allocator);
                if (s > 0) {
                    shiftLeftInPlace(v.limbs[0..v.size], &v.limbs[v.size], @intCast(s));
                    shiftLeftInPlace(r.limbs[0..r.size], &r.limbs[r.size], @intCast(s));
                }

                var j: u32 = m;
                while (true) : (j -= 1) {
                    const u_high: u32 = r.limbs[j + n];
                    const u_next: u32 = r.limbs[j + n - 1];
                    const qhat: u32 = estimateQuotientDigit(u_high, u_next, v.limbs[n - 1]);
                    o.limbs[j] = qhat;

                    if (mulSub(r.limbs[j .. j + n + 1], v.limbs[0..n], qhat)) {
                        o.limbs[j] -= 1;
                        addBack(r.limbs[j .. j + n], &r.limbs[j + n], v.limbs[0..n]);
                    }

                    if (j == 0) break;
                }

                o.size = m + 1;
                o.trimSize();

                if (s > 0) shiftRightInPlace(r.limbs[0..r.size], @intCast(s));

                o.positive = x.positive == y.positive;
            },
            .float => {
                // Check y is integer -> y == floor(y)
                // If not, use rational.add
            },
            .int => {},
            else => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .integer => {
                // Check x is integer -> x == floor(x)
            },
            else => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .integer => {},
            else => unreachable,
        },
        else => unreachable,
    }
}

fn shiftLeftInPlace(limbs: []u32, next: ?*u32, s: u5) void {
    if (s == 0 or limbs.len == 0) return;

    var carry: u32 = 0;
    for (limbs) |*limb| {
        const new_carry: u32 = limb.* >> @as(u5, @intCast(32 - @as(u6, @intCast(s))));
        limb.* = (limb.* << s) | carry;
        carry = new_carry;
    }

    if (carry != 0) {
        // Extend by one limb if necessary
        if (next) |n| {
            n.* = carry;
        }
    }
}

fn shiftRightInPlace(limbs: []u32, s: u5) void {
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

fn estimateQuotientDigit(u_high: u32, u_next: u32, v_high: u32) u32 {
    const B: u64 = 1 << 32;
    const num: u64 = (@as(u64, u_high) << 32) | u_next;
    var qhat: u64 = num / @as(u64, v_high);
    if (qhat > (B - 1)) qhat = B - 1;
    return @truncate(qhat);
}

fn mulSub(u_part: []u32, v: []const u32, qhat: u32) bool {
    var carry: u64 = 0;
    var borrow: u64 = 0;

    const n = v.len;
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const p: u64 = @as(u64, v[i]) * qhat + carry;
        const sub: u64 = @as(u64, u_part[i]) -% (p & 0xFFFFFFFF) -% borrow;
        u_part[i] = @truncate(sub);

        carry = p >> 32;
        borrow = if (sub >> 63 != 0) 1 else 0; // borrow if underflow
    }

    // subtract the final carry and borrow from the high limb
    const sub_high: u64 = @as(u64, u_part[n]) -% carry -% borrow;
    u_part[n] = @truncate(sub_high);

    return (sub_high >> 63) != 0; // true if borrow propagated out
}

fn addBack(u_part: []u32, next: *u32, v: []const u32) void {
    var carry: u64 = 0;
    const n = v.len;
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const sum: u64 = @as(u64, u_part[i]) + v[i] + carry;
        u_part[i] = @truncate(sum);
        carry = sum >> 32;
    }

    if (carry != 0) {
        next.* +%= @truncate(carry);
    }
}
