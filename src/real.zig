const std = @import("std");
const Rational = @import("rational.zig").Rational;

pub const Real = struct {
    rational: Rational,
    //irrational: []Irrational,

    pub const empty: Real = .{
        .rational = .empty,
    };
};
