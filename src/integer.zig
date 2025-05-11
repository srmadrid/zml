const std = @import("std");

pub const Integer = struct {
    limbs: []usize,
    metadata: usize,

    pub const empty: Integer = .{
        .limbs = &.{},
        .metadata = 0,
    };
};

pub const IntegerManaged = struct {
    limbs: []usize,
    metadata: usize,
    allocator: std.mem.Allocator,

    pub const empty: IntegerManaged = .{
        .limbs = &.{},
        .metadata = 0,
        .allocator = undefined,
    };
};
