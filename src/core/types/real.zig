const std = @import("std");
const RationalUnmanaged = @import("rational.zig").RationalUnmanaged;

pub const RealUnmanaged = struct {
    rational: RationalUnmanaged,
    //irrational: []Irrational,
};

pub const Real = struct {
    rational: RationalUnmanaged,
    //irrational: []Irrational,
    allocator: std.mem.Allocator,
};
