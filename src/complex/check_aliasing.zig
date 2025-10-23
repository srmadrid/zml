const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

pub fn check_aliasing(allocator: std.mem.Allocator, o: anytype, x: anytype) !@TypeOf(x) {
    var tx: @TypeOf(x) = undefined;
    if (comptime types.Scalar(types.Child(@TypeOf(x))) == rational.Rational) {
        if (comptime types.Scalar(@TypeOf(x)) == rational.Rational) { // Rational, Rational
            tx.re.num = if (o.re.num.limbs == x.re.num.limbs or o.re.den.limbs == x.re.num.limbs or
                o.im.num.limbs == x.re.num.limbs or o.im.den.limbs == x.re.num.limbs)
                try x.re.num.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.re.num;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.re.num.deinit(allocator);
            tx.re.den = if (o.re.num.limbs == x.re.den.limbs or o.re.den.limbs == x.re.den.limbs or
                o.im.num.limbs == x.re.den.limbs or o.im.den.limbs == x.re.den.limbs)
                try x.re.den.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.re.den;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.re.den.deinit(allocator);
            tx.im.num = if (o.re.num.limbs == x.im.num.limbs or o.re.den.limbs == x.im.num.limbs or
                o.im.num.limbs == x.im.num.limbs or o.im.den.limbs == x.im.num.limbs)
                try x.im.num.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.im.num;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.im.num.deinit(allocator);
            tx.im.den = if (o.re.num.limbs == x.im.den.limbs or o.re.den.limbs == x.im.den.limbs or
                o.im.num.limbs == x.im.den.limbs or o.im.den.limbs == x.im.den.limbs)
                try x.im.den.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.im.den;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
        } else { // Rational, Real
            tx.re.rational.num = if (o.re.num.limbs == x.re.rational.num.limbs or o.re.den.limbs == x.re.rational.num.limbs or
                o.im.num.limbs == x.re.rational.num.limbs or o.im.den.limbs == x.re.rational.num.limbs)
                try x.re.rational.num.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.re.rational.num;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.re.rational.num.deinit(allocator);
            tx.re.rational.den = if (o.re.num.limbs == x.re.rational.den.limbs or o.re.den.limbs == x.re.rational.den.limbs or
                o.im.num.limbs == x.re.rational.den.limbs or o.im.den.limbs == x.re.rational.den.limbs)
                try x.re.rational.den.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.re.rational.den;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.re.rational.den.deinit(allocator);
            std.debug.print("Also check for irrationals for tx.re.irrationals (when implemented)\n", .{});
            tx.im.rational.num = if (o.re.num.limbs == x.im.rational.num.limbs or o.re.den.limbs == x.im.rational.num.limbs or
                o.im.num.limbs == x.im.rational.num.limbs or o.im.den.limbs == x.im.rational.num.limbs)
                try x.im.rational.num.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.im.rational.num;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.im.rational.num.deinit(allocator);
            tx.im.rational.den = if (o.re.num.limbs == x.im.rational.den.limbs or o.re.den.limbs == x.im.rational.den.limbs or
                o.im.num.limbs == x.im.rational.den.limbs or o.im.den.limbs == x.im.rational.den.limbs)
                try x.im.rational.den.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.im.rational.den;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            std.debug.print("Also check for irrationals for tx.im.irrationals (when implemented)\n", .{});
        }
    } else {
        if (comptime types.Scalar(@TypeOf(x)) == rational.Rational) { // Real, Rational
            tx.re.num = if (o.re.rational.num.limbs == x.re.num.limbs or o.re.rational.den.limbs == x.re.num.limbs or
                o.im.rational.num.limbs == x.re.num.limbs or o.im.rational.den.limbs == x.re.num.limbs)
                try x.re.num.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.re.num;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.re.num.deinit(allocator);
            tx.re.den = if (o.re.rational.num.limbs == x.re.den.limbs or o.re.rational.den.limbs == x.re.den.limbs or
                o.im.rational.num.limbs == x.re.den.limbs or o.im.rational.den.limbs == x.re.den.limbs)
                try x.re.den.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.re.den;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.re.den.deinit(allocator);
            tx.im.num = if (o.re.rational.num.limbs == x.im.num.limbs or o.re.rational.den.limbs == x.im.num.limbs or
                o.im.rational.num.limbs == x.im.num.limbs or o.im.rational.den.limbs == x.im.num.limbs)
                try x.im.num.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.im.num;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.im.num.deinit(allocator);
            tx.im.den = if (o.re.rational.num.limbs == x.im.den.limbs or o.re.rational.den.limbs == x.im.den.limbs or
                o.im.rational.num.limbs == x.im.den.limbs or o.im.rational.den.limbs == x.im.den.limbs)
                try x.im.den.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.im.den;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
        } else { // Real, Real
            tx.re.rational.num = if (o.re.rational.num.limbs == x.re.rational.num.limbs or o.re.rational.den.limbs == x.re.rational.num.limbs or
                o.im.rational.num.limbs == x.re.rational.num.limbs or o.im.rational.den.limbs == x.re.rational.num.limbs)
                try x.re.rational.num.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.re.rational.num;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.re.rational.num.deinit(allocator);
            tx.re.rational.den = if (o.re.rational.num.limbs == x.re.rational.den.limbs or o.re.rational.den.limbs == x.re.rational.den.limbs or
                o.im.rational.num.limbs == x.re.rational.den.limbs or o.im.rational.den.limbs == x.re.rational.den.limbs)
                try x.re.rational.den.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.re.rational.den;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.re.rational.den.deinit(allocator);
            std.debug.print("Also check for irrationals for tx.re.irrationals (when implemented)\n", .{});
            tx.im.rational.num = if (o.re.rational.num.limbs == x.im.rational.num.limbs or o.re.rational.den.limbs == x.im.rational.num.limbs or
                o.im.rational.num.limbs == x.im.rational.num.limbs or o.im.rational.den.limbs == x.im.rational.num.limbs)
                try x.im.rational.num.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.im.rational.num;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            errdefer tx.im.rational.num.deinit(allocator);
            tx.im.rational.den = if (o.re.rational.num.limbs == x.im.rational.den.limbs or o.re.rational.den.limbs == x.im.rational.den.limbs or
                o.im.rational.num.limbs == x.im.rational.den.limbs or o.im.rational.den.limbs == x.im.rational.den.limbs)
                try x.im.rational.den.copy(allocator)
            else blk: {
                var tmp: integer.Integer = x.im.rational.den;
                tmp.flags.owns_data = false;
                break :blk tmp;
            };
            std.debug.print("Also check for irrationals for tx.im.irrationals (when implemented)\n", .{});
        }
    }

    tx.flags = .{ .owns_data = true, .writable = false };

    return tx;
}
