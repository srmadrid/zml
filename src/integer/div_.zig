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
                if (y.size == 0) return integer.Error.ZeroDivision;

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
                    o.trimSize();
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
                    shiftLeftInPlace(v.limbs[0..v.size], null, @intCast(s));
                    shiftLeftInPlace(r.limbs[0..r.size], &r.limbs[r.size], @intCast(s));
                }

                try o.reserve(allocator, m + 1);

                var j: u32 = m;
                while (true) : (j -= 1) {
                    const qhat: u32 = estimate_qhat(
                        r.limbs[j + n],
                        r.limbs[j + n - 1],
                        r.limbs[j + n - 2],
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
                o.trimSize();

                if (s > 0) shiftRightInPlace(r.limbs[0..n], @intCast(s));
                if (r.limbs[r.size - 1] == 0 and r.size > 1) r.size -= 1;

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

pub fn shiftLeftInPlace(limbs: []u32, next: ?*u32, s: u5) void {
    if (s == 0 or limbs.len == 0) return;

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
    }
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

        const res = @subWithOverflow(@as(u64, u_part[i]), (prod & 0xFFFFFFFF) + borrow);
        u_part[i] = @truncate(res[0]);
        borrow = @intCast(res[1]);
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
