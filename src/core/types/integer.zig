const std = @import("std");

pub const IntegerUnmanaged = struct {
    limbs: []usize,
    metadata: usize,

    pub const empty: IntegerUnmanaged = .{
        .limbs = &.{},
        .metadata = 0,
    };
};

pub const Integer = struct {
    limbs: []usize,
    metadata: usize,
    allocator: std.mem.Allocator,

    pub const empty: Integer = .{
        .limbs = &.{},
        .metadata = 0,
        .allocator = undefined,
    };
};
