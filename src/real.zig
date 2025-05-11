const std = @import("std");
const Rational = @import("rational.zig").Rational;

pub const Real = struct {
    rational: Rational,
    //irrational: []Irrational,

    pub const empty: Real = .{
        .rational = .empty,
    };
};

pub const RealManaged = struct {
    rational: Rational,
    //irrational: []Irrational,
    allocator: std.mem.Allocator,

    pub const empty: RealManaged = .{
        .rational = .empty,
        .allocator = undefined,
    };
};
