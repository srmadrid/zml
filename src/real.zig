const std = @import("std");
const Rational = @import("rational.zig").Rational;

pub const Real = struct {
    rational: Rational,
    //irrationals: []Irrational, Maybe istead hold something like struct { irrational: Irrational, multiplicity: i32 }, with for instance, .{pi, -2} means \pi^{-2}

    pub const empty: Real = .{
        .rational = .empty,
    };
};
