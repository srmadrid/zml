const std = @import("std");

const zsl = @import("zsl");

const tzsl = @import("zsl.zig");

pub fn Static(N: type) type {
    return struct {
        pub const instantiate = true;
        pub const Numeric = N;
    };
}

pub fn printVector(desc: []const u8, v: anytype) void {
    std.debug.print("\nVector {s}:\n", .{desc});

    const v_len = if (comptime zsl.meta.isStaticVector(@TypeOf(v))) @TypeOf(v).len else v.len;

    var i: usize = 0;
    while (i < v_len) : (i += 1) {
        if (comptime zsl.meta.isComplex(zsl.meta.Numeric(@TypeOf(v)))) {
            std.debug.print("{d} + {d}i\n", .{ (v.getAssumeInBounds(i)).re, (v.getAssumeInBounds(i)).im });
        } else {
            std.debug.print("{d}\n", .{v.getAssumeInBounds(i)});
        }
    }
    std.debug.print("\n", .{});
}

pub fn randomVector(comptime V: type, allocator: std.mem.Allocator, rand: std.Random, comptime len: usize, inc: isize, nnz: usize) !V {
    switch (comptime zsl.meta.vectorType(V)) {
        .static => {
            var result: V = .init;

            inline for (0..len) |i| {
                result.setAssumeInBounds(i, tzsl.randomNumber(zsl.meta.Numeric(V), rand));
            }

            return result;
        },
        .dense => {
            var result: V = try .init(allocator, len * zsl.numeric.cast(usize, zsl.int.abs(inc)));
            result.flags.noconj = rand.boolean();

            var i: usize = 0;
            while (i < result.len) : (i += 1) {
                result.setAssumeInBounds(i, tzsl.randomNumber(zsl.meta.Numeric(V), rand));
            }

            result.len = len;
            result.inc = inc;

            return result;
        },
        .sparse => {
            var result: V = try .init(allocator, len, nnz);
            errdefer result.deinit(allocator);
            result.flags.noconj = rand.boolean();

            // generate random indices
            var count: usize = 0;
            while (count < nnz) : (count += 1) {
                const i = rand.intRangeAtMost(usize, 0, len - 1);
                try result.set(allocator, i, tzsl.randomNumber(zsl.meta.Numeric(V), rand));
            }

            return result;
        },
        else => unreachable,
    }
}

pub fn correctApply2(comptime O: type, allocator: std.mem.Allocator, len: usize, u: anytype, v: anytype, op_: anytype) !zsl.vector.Dense(O) {
    const result: zsl.vector.Dense(O) = try .init(allocator, len);

    var i: usize = 0;
    while (i < result.len) : (i += 1) {
        op_(
            &result.data[result._index(i)],
            if (comptime zsl.meta.isVector(@TypeOf(u))) u.get(i) catch unreachable else u,
            if (comptime zsl.meta.isVector(@TypeOf(v))) v.get(i) catch unreachable else v,
        );
    }

    return result;
}

pub fn areEql(u: anytype, v: anytype) !void {
    const u_len = switch (comptime zsl.meta.vectorType(@TypeOf(u))) {
        .static => @TypeOf(u).len,
        .dense, .sparse => u.len,
        else => unreachable,
    };

    if (u_len != v.len)
        return error.NotEqual;

    var all_eql = true;

    var i: usize = 0;
    while (i < v.len) : (i += 1) {
        all_eql = all_eql and zsl.numeric.eq(u.get(i) catch unreachable, v.get(i) catch unreachable);
    }

    if (!all_eql)
        return error.NotEqual;
}

test {
    const test_apply2 = false;

    if (test_apply2) {
        _ = @import("vector/apply2.zig");
        _ = @import("vector/apply2Alloc.zig");
        _ = @import("vector/apply2Into.zig");
    }
}
