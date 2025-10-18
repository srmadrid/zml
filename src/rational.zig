const std = @import("std");

const types = @import("types.zig");
const ops = @import("ops.zig");
const int = @import("int.zig");
const integer = @import("integer.zig");
const Integer = integer.Integer;

pub const Rational = struct {
    num: Integer,
    den: Integer,
    flags: Flags,

    pub const empty: Rational = .{
        .num = .empty,
        .den = .empty,
    };

    pub fn init(allocator: std.mem.Allocator, numsize: u32, densize: u32) !Rational {
        var num: Integer = try Integer.init(allocator, numsize);
        errdefer num.deinit(allocator);

        var den: Integer = try Integer.init(allocator, densize);
        errdefer den.deinit(allocator);

        return .{
            .num = num,
            .den = den,
            .flags = .{ .owns_data = true, .writable = true },
        };
    }

    pub fn initSet(allocator: std.mem.Allocator, numerator: anytype, denominator: anytype) !Rational {
        // Edit to properly handle floats, cfloats, complexes and other rationals correctly
        if (ops.eq(denominator, 0, .{}) catch unreachable)
            return Error.ZeroDenominator;

        var num: Integer = try Integer.initSet(allocator, numerator);
        errdefer num.deinit(allocator);

        var den: Integer = try Integer.initSet(allocator, denominator);
        errdefer den.deinit(allocator);

        var r: Rational = .{
            .num = num,
            .den = den,
            .flags = .{ .owns_data = true, .writable = true },
        };

        try r.reduce(allocator);

        return r;
    }

    pub fn deinit(self: *Rational, allocator: std.mem.Allocator) void {
        if (self.flags.owns_data) {
            self.num.deinit(allocator);
            self.den.deinit(allocator);
        }

        self.* = undefined;
    }

    pub fn reduce(self: *Rational, allocator: std.mem.Allocator) !void {
        if (!self.flags.writable)
            return Error.NotWritable;

        var g: Integer = try integer.gcd(allocator, self.num, self.den);
        defer g.deinit(allocator);

        try integer.div_(allocator, &self.num, self.num, g);
        try integer.div_(allocator, &self.den, self.den, g);
    }

    pub fn copy(self: *const Rational, allocator: std.mem.Allocator) !Rational {
        var num: Integer = try self.num.copy(allocator);
        errdefer num.deinit(allocator);
        const den: Integer = try self.den.copy(allocator);

        return .{
            .num = num,
            .den = den,
            .flags = .{ .owns_data = true, .writable = true },
        };
    }
};

// Arithmetic operations
// pub const add = @import("rational/add.zig").add;
pub const add_ = @import("rational/add_.zig").add_;
// pub const sub = @import("rational/sub.zig").sub;
pub const sub_ = @import("rational/sub_.zig").sub_;
// pub const mul = @import("rational/mul.zig").mul;
pub const mul_ = @import("rational/mul_.zig").mul_;
// pub const div = @import("rational/div.zig").div;
pub const div_ = @import("rational/div_.zig").div_;

// Comparison operations
// pub const cmp = @import("rational/cmp.zig").cmp;
// pub const eq = @import("rational/eq.zig").eq;
// pub const ne = @import("rational/ne.zig").ne;
// pub const lt = @import("rational/lt.zig").lt;
// pub const le = @import("rational/le.zig").le;
// pub const gt = @import("rational/gt.zig").gt;
// pub const ge = @import("rational/ge.zig").ge;

// Basic operations
pub const abs = @import("rational/abs.zig").abs;
pub const neg = @import("rational/neg.zig").neg;

pub const Error = error{
    ZeroDenominator,
    ZeroDivision,
    NotWritable,
};

pub const Flags = packed struct {
    owns_data: bool = true,
    writable: bool = true,
};
