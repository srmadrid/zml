const std = @import("std");

const types = @import("types.zig");
const ops = @import("ops.zig");
const integer = @import("integer.zig");
const Integer = integer.Integer;

pub const Rational = struct {
    num: Integer,
    den: Integer,

    pub const empty: Rational = .{
        .num = .empty,
        .den = .empty,
    };

    pub fn initSet(allocator: std.mem.Allocator, numerator: anytype, denominator: anytype) !Rational {
        if (ops.eq(denominator, 0, .{}) catch unreachable)
            return Error.ZeroDenominator;

        var num: Integer = try Integer.initSet(allocator, numerator);
        errdefer num.deinit(allocator);

        var den: Integer = try Integer.initSet(allocator, denominator);
        errdefer den.deinit(allocator);

        var r: Rational = .{
            .num = num,
            .den = den,
        };

        try r.reduce(allocator);

        return r;
    }

    pub fn deinit(self: *Rational, allocator: std.mem.Allocator) void {
        self.num.deinit(allocator);
        self.den.deinit(allocator);
    }

    pub fn reduce(self: *Rational, allocator: std.mem.Allocator) !void {
        std.debug.print("Reducing rational\n", .{});
        var g: Integer = try integer.gcd(allocator, self.num, self.den);
        defer g.deinit(allocator);
        std.debug.print("Extracted GCD\n", .{});

        try integer.div_(allocator, &self.num, self.num, g);
        std.debug.print("Reduced numerator\n", .{});
        try integer.div_(allocator, &self.den, self.den, g);
        std.debug.print("Reduced denominator\n", .{});
    }
};

pub const Error = error{
    ZeroDenominator,
};
